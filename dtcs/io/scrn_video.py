"""TODO(Andrew)"""

import copy
import logging
import math
import os

_logger = logging.getLogger(__name__)

import matplotlib.pyplot as plt
import numpy as np

pygame = False
try:
    import pygame
    from dtcs.sim.surface_crn.surface_crns.views.grid_display import (
        HexGridDisplay, ParallelEmulatedSquareGridDisplay, SquareGridDisplay)
    from dtcs.sim.surface_crn.surface_crns.views.coord_grid_display import \
        CoordGridDisplay
    from dtcs.sim.surface_crn.hex_grid_with_intersect import HexGridPlusIntersectionDisplay
    from dtcs.sim.surface_crn.surface_crns.views.legend_display import LegendDisplay
    from dtcs.sim.surface_crn.surface_crns.views.text_display import TextDisplay
    from dtcs.sim.surface_crn.surface_crns.views.time_display import TimeDisplay

    from dtcs.sim.surface_crn.surface_crns.simulators.queue_simulator import QueueSimulator
    from dtcs.sim.surface_crn.surface_crns.pygbutton import *
except ModuleNotFoundError:
    _logger.info('Didn\'t load module pygame')

try: import cv2
except ModuleNotFoundError: _logger.info('Didn\'t load module cv2')

MOVIE_SUBDIRECTORY = 'movies'
FRAME_SUBDIRECTORY = 'frames'

PROFILE = False
WHITE = (255, 255, 255)
BLACK = (0, 0, 0)
# time_font = pygame.font.SysFont('monospace', 24)
# TODO: unify font with other figures
# time_font = time_font = pygame.font.SysFont(pygame.font.get_default_font(), 24)
# import matplotlib
# time_font = time_font = pygame.font.SysFont(matplotlib.rcParams['font.family'], 24)
TEXT_X_BUFFER = 10
TEXT_Y_BUFFER = 5
if pygame:
    time_font = pygame.font.SysFont('arialttf', 24)
    TEXT_HEIGHT = time_font.get_linesize() + 2 * TEXT_Y_BUFFER
button_width = 60
button_height = 30
button_buffer = 5
MIN_GRID_WIDTH = 6 * button_width + 10 * button_buffer

def _get_surface_display(surface, init_surface, opts, display_class=None):
    if opts.grid_type == 'parallel_emulated':
        return ParallelEmulatedSquareGridDisplay(
            grid=init_surface,
            colormap=opts.COLORMAP,
            emulation_colormap=opts.emulation_colormap,
            horizontal_buffer=opts.horizontal_buffer,
            vertical_buffer=opts.vertical_buffer,
            cell_height=opts.cell_height,
            cell_width=opts.cell_width,
            representative_cell_x=opts.representative_cell_x,
            representative_cell_y=opts.representative_cell_y,
            min_x=MIN_GRID_WIDTH,
            min_y=0,
            pixels_per_node=opts.pixels_per_node,
            display_text=opts.display_text
        )
    elif opts.grid_type == 'standard':
        if display_class:
            DisplayClass = display_class
        elif surface.structure == 'voronoi':
            DisplayClass = CoordGridDisplay
        elif surface.structure == 'hexagon':
            DisplayClass = HexGridPlusIntersectionDisplay
        elif opts.grid_type == "standard" and opts.surface_geometry == "square":
            DisplayClass = SquareGridDisplay
        elif opts.grid_type == "standard" and opts.surface_geometry == "hex":
            DisplayClass = HexGridDisplay
        else:
            assert False, 'No display class but grid type is standard'
        return DisplayClass(grid=init_surface,
                            colormap=opts.COLORMAP,
                            min_x=MIN_GRID_WIDTH,
                            min_y=0,
                            pixels_per_node=opts.pixels_per_node,
                            display_text=opts.display_text)
    else:
        raise Exception("Unrecognized grid type '" + opts.grid_type + "'")

def _make_imgbuf_from_plot(
        plot_func,
        args=None,
        kwargs=None
):
    args = args or []
    kwargs = kwargs or {}

    fig = plot_func(*args, **kwargs)

    fig.canvas.draw()
    plt.close()

    image_rgb256 = np.frombuffer(fig.canvas.tostring_rgb(), dtype=np.uint8)
    image_rgb256 = image_rgb256.reshape(fig.canvas.get_width_height()[::-1] + (3,))

    return image_rgb256, (image_rgb256.shape[1], image_rgb256.shape[0])

def _make_video(frames, display_size, fps):

    try:
        out = cv2.VideoWriter(
            filename='output.mp4',
            fourcc=cv2.VideoWriter_fourcc(*'LMP4'),
            fps=fps,
            frameSize=display_size,
        )

        for filename in frames:
            img = cv2.imread(filename)
            out.write(img)
    finally:
        out.release()

    return os.path.abspath('output.mp4')


def _make_scrn_frames(
        scrn,
        run=0,

        plot=False,
        plot_kwargs=None,

        frame_times=(0, ),
        frames_dir='.',
):
    plot_kwargs = plot_kwargs or {}
    frame_time_queue = list(frame_times)

    # --- Prepare to save output ------------------------------------
    for directory in [frames_dir]:
        if not os.path.isdir(directory):
            os.mkdir(directory)

    # --- Prepare movie settings --------------------------------------
    # TODO(Andrew): move to make_movie
    # frame_time_queue = list(np.arange(0, scrn.time_max, 1 / frames_per_timestep))

    # --- Retreive run information -----------------------------------------
    # Retreive information about that run
    event_history = scrn._event_histories[run]
    last_event = event_history[-1]
    surface = copy.deepcopy(scrn._init_surfaces[run])

    # Make a dummy simulation
    simulator = QueueSimulator(
        surface=surface,
        rxns=scrn.rsys,  # TODO(Andrew) Can I remove this? What's it for?
    )

    # def count_species(): return dict(simulator.surface.species_count())
    # times = [0]
    # species_counts = [count_species()]

    # --- Set up pygame ------------------------------------------------
    grid_display = _get_surface_display(
        scrn.crn.surface,
        surface,
        scrn.crn,
        display_class=None
    )
    scrn.crn.debug = True

    # #############################################################################
    # --- Create Displays for Legend, Title, etc. ---
    legend_display = LegendDisplay(colormap=scrn.crn.COLORMAP)
    # Width only requires legend and grid sizes to calculate
    display_width = grid_display.display_width + legend_display.display_width
    # Width used to calculate time label and button placements
    time_display = TimeDisplay(display_width)
    # TODO(Andrew): Add back the title

    # Display for the additional title
    title_display = TextDisplay(display_width, text="Surface CRN Trajectory")
    show_title = bool(plot)

    display_height = max(
        legend_display.display_height + 2 * legend_display.VERTICAL_BUFFER,
        grid_display.display_height + 1
    )
    display_height += time_display.display_height
    if show_title:
        display_height += title_display.display_height

    # Adjust display for additional plot
    if callable(plot):
        plot_buf, (plot_width, plot_height) = _make_imgbuf_from_plot(
            plot_func=plot,
            kwargs=dict(
                scts=scrn,
                run=run,
                time=0,
                **plot_kwargs,
            )
        )

        plot_xpos = display_width
        plot_ypos = 0

        display_width = plot_xpos + plot_width
        display_height = max(display_height, plot_ypos + plot_height)

    # --- Initialize PyGame Canvas --------------------------------------------
    display_surface = pygame.Surface(
        (display_width, display_height),
        flags=0,
        depth=32,
    )
    # display_surface = pygame.display.set_mode((display_width, display_height), 0, 32)
    simulator.display_surface = display_surface

    display_surface.fill(WHITE)

    # --- Blit Static Features (e.g. title) onto the Canvas -------------------
    # these do the rendering on the surface
    if show_title:
        title_display.render(display_surface, x_pos=0,
                             y_pos=0)
        time_display.render(display_surface, x_pos=0,
                            y_pos=title_display.y_pos +
                                  title_display.display_height)
    else:
        time_display.render(display_surface, x_pos=0, y_pos=0)
    legend_display.render(display_surface, x_pos=0,
                          y_pos=time_display.y_pos +
                                time_display.display_height)

    # legend_display.display_height is intended to be a very small number.
    grid_display.render(
        display_surface,
        x_pos=legend_display.display_width,
        y_pos=time_display.y_pos + time_display.display_height,
        width=time_display.display_width - legend_display.display_width,
        height=legend_display.display_height * 4
    )
    # ##########################################################################

    # --- Iterate through the simulation and make frames ----------------------
    frame_filepaths = []
    for reaction in event_history:
        # If we have made all the frames, we're done
        if not frame_time_queue:
            break

        # If the current reaction time is after the next frame, make a frame
        #  and remove that frame time from the queue.
        # Also make a frame if this is the last event, and there's still frames
        #  remaining.
        if reaction.time >= frame_time_queue[0] or reaction == last_event:
            frame_time = frame_time_queue.pop(0)

            # Update the time display
            time_display.time = frame_time
            time_display.render(
                display_surface,
                x_pos=time_display.x_pos,
                y_pos=time_display.y_pos
            )

            # Render the plotting function
            if plot:
                plot_buf, _ = _make_imgbuf_from_plot(
                    plot_func=plot,
                    kwargs=dict(
                        scts=scrn,
                        run=run,
                        time=frame_time,
                        **plot_kwargs,
                    )
                )
                plot_image = pygame.image.frombuffer(
                    plot_buf,
                    (plot_width, plot_height),
                    'RGB'
                )
                display_surface.blit(
                    plot_image,
                    (plot_xpos, plot_ypos),
                )

            # Re-render the whole frame
            #  It seems like this would take more time, but it actually might not
            #  or it might. It's easier to code, lol
            grid_display.render(
                display_surface,
                x_pos=legend_display.display_width,
                y_pos=time_display.y_pos + time_display.display_height,
                width=time_display.display_width - legend_display.display_width,
                height=legend_display.display_height * 4
            )

            # Do all the updating before outputting
            # pygame.display.update()
            # pygame.display.flip()

            # Save the frame
            filepath = os.path.join(frames_dir, f'frame_t{frame_time}.png')
            frame_filepaths.append(filepath)
            pygame.image.save(display_surface, filepath)
            _logger.debug(f'Worte frame {filepath}')
            # TODO(Andrew): Print what frames are being made as debug

        # Manually apply the reaction
        for index in range(len(reaction.participants)):
            cell = surface.get_node_by_id(reaction.participants[index].node_id)
            state = reaction.rule.outputs[index]
            cell.state = state

    return frame_filepaths, display_width, display_height


def make_scrn_image(
        scrn,
        run=0,
        time=-1,

        plot=False,
        plot_kwargs=None,

        output_dir='.',
):
    # --- Prepare to save output ------------------------------------
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    # --- Sanitize frame time ---------------------------------------
    if time < 0:
        time = scrn.time_max

    # --- Make the image -------------------------------------------
    frame_filepaths, frame_width, frame_height = _make_scrn_frames(
        scrn=scrn,
        run=run,

        plot=plot,
        plot_kwargs=plot_kwargs,

        frame_times=[time],
        frames_dir=output_dir,
    )

    return frame_filepaths[0]


def make_scrn_video(
        scrn,
        run=0,
        plot=False,
        plot_kwargs=None,

        frames_per_timestep=10,
        output_dir='output',
):
    # --- Prepare to save output ----------------------------------------------
    movie_dir = output_dir
    frames_dir = os.path.join(output_dir)
    for directory in [output_dir, movie_dir, frames_dir]:
        if not os.path.isdir(directory):
            os.mkdir(directory)

    # --- Prepare frame times -------------------------------------------------
    time_precision = round(math.log(frames_per_timestep, 10) + 1)
    frame_times = list(np.arange(0, scrn.time_max, 1 / frames_per_timestep))
    frame_times = list(map(
        lambda t: round(t, time_precision),
        frame_times,
    ))

    frame_filepaths, frame_width, frame_height = _make_scrn_frames(
        scrn=scrn,
        run=run,

        plot=plot,
        plot_kwargs=plot_kwargs,

        frame_times=frame_times,
        frames_dir=frames_dir,
    )

    return _make_video(
        frame_filepaths,
        (frame_width, frame_height),
        fps=2,
    )


# from dtcs.sim.surface_crn.surface_crns.SurfaceCRNQueueSimulator import display_next_event

# def _update_display(
#     opts,
#     simulator,
#     time,
#     grid_display,
#     time_display=None,
#     title_display=None,
# ):
#     grid_display.re_render()
#
# def _save_frame():
#     pass
#
# import numpy as np
