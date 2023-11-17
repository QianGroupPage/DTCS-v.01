"""TODO(Andrew)"""

from typing import Optional
import copy
import logging
import math
import os

import matplotlib.pyplot as plt
import numpy as np

from dtcs.common import util


_logger = logging.getLogger(__name__)


# --- Optional Imports --------------------------------------------------------
try: import pygame
except ModuleNotFoundError: _logger.info('Didn\'t load module pygame')
try: import cv2
except ModuleNotFoundError: _logger.info('Didn\'t load module cv2')
try: from PIL import Image
except ModuleNotFoundError: _logger.info('Didn\'t load module Pillow')

if util.feature_loaded('scrn-image'):
    from dtcs.sim.surface_crn.surface_crns.views.grid_display import (
        HexGridDisplay, ParallelEmulatedSquareGridDisplay, SquareGridDisplay)
    from dtcs.sim.surface_crn.surface_crns.views.coord_grid_display import \
        CoordGridDisplay
    from dtcs.sim.surface_crn.hex_grid_with_intersect import HexGridPlusIntersectionDisplay
    from dtcs.sim.surface_crn.surface_crns.simulators.queue_simulator import QueueSimulator

MOVIE_SUBDIRECTORY = 'movies'
FRAME_SUBDIRECTORY = 'frames'
MIN_GRID_WIDTH = 410


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


def _make_video(frames, fps, out_file='output.mp4'):
    frame_dims = Image.open(frames[0]).size

    out = cv2.VideoWriter(
        filename=out_file,
        fourcc=cv2.VideoWriter_fourcc(*'LMP4'),
        fps=fps,
        frameSize=frame_dims,
    )

    for filename in frames:
        img = cv2.imread(filename)
        out.write(img)
    out.release()

    return os.path.abspath(out_file)


def make_scrn_images(
        scts,
        frame_times,
        run=0,
        dpi=100,
):
    frame_time_queue = sorted(frame_times)

    # --- Retreive run information -----------------------------------------
    # Retreive information about that run
    event_history = scts._event_histories[run]
    last_event = event_history[-1]
    surface = copy.deepcopy(scts._init_surfaces[run])

    # Make a dummy simulation
    simulator = QueueSimulator(
        surface=surface,
        rxns=scts.rsys,  # TODO(Andrew) Can I remove this? What's it for?
    )

    # --- Get the display of the surface --------------------------------------
    surface_display = _get_surface_display(
        scts.crn.surface,
        surface,
        scts.crn,
        display_class=None
    )
    surface_display.fig_dpi = dpi
    surface_display.recalculate_display_sizes()
    display_width = surface_display.display_width
    display_height = surface_display.display_height

    # --- Initialize PyGame Canvas --------------------------------------------
    canvas = pygame.Surface(
        (display_width, display_height),
        flags=0,
        depth=32,
    )
    # canvas = pygame.display.set_mode((display_width, display_height), 0, 32)
    simulator.display_surface = canvas
    canvas.fill((255, 255, 255))

    # --- Iterate through the simulation and make frames ----------------------
    images = {}
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

            # Render the whole frame
            surface_display.render(
                canvas,
                x_pos=0,
                y_pos=0,
                width=display_width,
                height=display_height,
            )

            # Save the frame
            img_str = pygame.image.tostring(canvas, 'RGBA')
            img_pil = Image.frombytes(
                mode='RGBA',
                size=(display_width, display_height),
                data=img_str,
            )
            images[frame_time] = np.array(img_pil)
            _logger.debug(f'Worte frame t={frame_time}')

        # Manually apply the reaction
        for index in range(len(reaction.participants)):
            cell = surface.get_node_by_id(reaction.participants[index].node_id)
            state = reaction.rule.outputs[index]
            cell.state = state

    return images, display_width, display_height


def make_scrn_video(
    scts,
    plot_func,
    run: int = 0,

    frames_dir: str = 'frames',  # TODO(Andrew) option to delete when done?
    output_fname: str = 'output',
    frames_per_timestep=10,
    frames_per_second=2,
    surface_img_dpi=200,

    **plot_kwargs
):
    # --- Prepare to save output ----------------------------------------------
    for directory in [frames_dir]:
        if not os.path.isdir(directory):
            os.mkdir(directory)

    # --- Prepare frame times -------------------------------------------------
    time_precision = round(math.log(frames_per_timestep, 10) + 1)
    frame_times = list(np.arange(0, scts.time_max, 1 / frames_per_timestep))
    frame_times = sorted(map(
        lambda t: round(t, time_precision),
        frame_times,
    ))

    # --- Make all the images at once -----------------------------------------
    images, width, height = make_scrn_images(
        scts=scts,
        run=run,
        frame_times=frame_times,
        dpi=surface_img_dpi,
    )

    # --- Plot each image and save to file ------------------------------------
    frame_filepaths = []
    for time, image in images.items():
        fig, axes = plot_func(
            surf_img=image,
            scts=scts,
            run=run,
            time=time,
            **plot_kwargs,
        )

        filepath = os.path.join(frames_dir, f'frame_t{time}.png')
        fig.savefig(filepath)
        plt.close()
        frame_filepaths.append(filepath)
        _logger.debug(f'Worte frame {filepath}')

    return _make_video(
        frame_filepaths,
        fps=frames_per_second,
        out_file=output_fname + '.mp4',
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
