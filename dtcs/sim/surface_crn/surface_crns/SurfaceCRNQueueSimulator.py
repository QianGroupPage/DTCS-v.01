"""
Simulates a surface chemical reaction network (CRN) on a 2D lattice.

Usage: python surface_CRN_simulator.py -m <manifest>

Simulates the stochastic behavior of a surface CRN on a 2D grid lattice. Only
implements unimolecular and bimolecular rules. Uses a Gillespie-Gibson-Bruck-
like algorithm for computing next reactions with a priority queue.
"""

from __future__ import print_function

try:
    import surface_crns
except ImportError:
    import sys

    sys.path.append("/")
import cProfile
import math
import optparse
import os
import subprocess as sp
import sys
from time import process_time

import dtcs.sim.surface_crn.surface_crns.readers as readers
from dtcs.common import ipython_visuals
from dtcs.sim.surface_crn.results import Results
from dtcs.sim.surface_crn.surface_crns.options.option_processor import \
    SurfaceCRNOptionParser
from dtcs.sim.surface_crn.surface_crns.pygbutton import *
from dtcs.sim.surface_crn.surface_crns.simulators.event_history import \
    EventHistory
from dtcs.sim.surface_crn.surface_crns.simulators.queue_simulator import \
    QueueSimulator
from dtcs.sim.surface_crn.surface_crns.simulators.synchronous_simulator import \
    SynchronousSimulator
from dtcs.sim.surface_crn.surface_crns.views.grid_display import (
    HexGridDisplay, ParallelEmulatedSquareGridDisplay, SquareGridDisplay)
from dtcs.sim.surface_crn.surface_crns.views.legend_display import LegendDisplay
from dtcs.sim.surface_crn.surface_crns.views.text_display import TextDisplay
from dtcs.sim.surface_crn.surface_crns.views.time_display import TimeDisplay

# TODO: study whether commenting these out would cause any unintended effects.
# These prevent pygame from opening when the entire library is loaded or when
# individual simulation is run.
# pygame.display.init()
# pygame.font.init()

#############
# CONSTANTS #
#############
PROFILE = False
WHITE = (255, 255, 255)
BLACK = (0, 0, 0)
# time_font = pygame.font.SysFont('monospace', 24)
# TODO: unify font with other figures
# time_font = time_font = pygame.font.SysFont(pygame.font.get_default_font(), 24)
# import matplotlib
# time_font = time_font = pygame.font.SysFont(matplotlib.rcParams['font.family'], 24)
time_font = time_font = pygame.font.SysFont('arialttf', 24)
TEXT_X_BUFFER = 10
TEXT_Y_BUFFER = 5
TEXT_HEIGHT = time_font.get_linesize() + 2 * TEXT_Y_BUFFER
button_width = 60
button_height = 30
button_buffer = 5
MIN_GRID_WIDTH = 6 * button_width + 10 * button_buffer

###########################
# MOVIE CAPTURE CONSTANTS #
###########################
MOVIE_SUBDIRECTORY = "movies"
DEBUG_SUBDIRECTORY = "debug"
FRAME_SUBDIRECTORY = "frames"
CUTOFF_TIME = 600000000  # Cut off simulation at 10 minutes
CUTOFF_SIZE = 10000 * 500000000000  # Cut off simulation at roughly 1000 frames for
# a typical image size.

#############
# VARIABLES #
#############
time = 0
simulation = None


################
# MAIN PROGRAM #
################

def main():
    available_options = optparse.OptionParser()
    available_options.add_option("-m", '--manifest', action="store",
                                 type='string', dest='manifest_filename')
    (command_line_options, args) = available_options.parse_args(sys.argv)
    manifest_filename = command_line_options.manifest_filename
    if not manifest_filename:
        raise Exception("Manifest file required (use the flag -m <filename>, " +
                        "where <filename> is the name of your manifest file)")
    simulate_surface_crn(manifest_filename)


def simulate_surface_crn(manifest_filename,
                         group_selection_seed,
                         display_class=None,
                         init_state=None,
                         rxns=None,
                         spectra_in_video=True,
                         running_average=10,
                         spectra_max_conc=-1,
                         ir_intensities=None):
    """
    Runs a simulation, and displays it in a GUI window OR saves all frames
    as PNG images.

    Normal operation is to read all options from a manifest file, given by
    manifest_filename. If you want to use a custom surface geometry (anything
    other than a square or hex grid), you'll need to supply your own initial
    state (whatever your display object will use, but must be an iterable and
    contain surface_crns.base.Node objects) and a class that can display your
    state (should subclass surface_crns.views.grid_display.SurfaceDisplay),
    which you should pass as "init_state" and "DisplayClass", respectively.
    """

    ################################
    # READ MANIFEST AND INITIALIZE #
    ################################
    # Parse the manifest
    # print("Reading information from manifest file " + manifest_filename + "...",
    #       end="")
    # TODO(Andrew): Replace this with CRNSpec
    manifest_options = \
        readers.manifest_readers.read_manifest(manifest_filename)
    opts = SurfaceCRNOptionParser(manifest_options)

    # print(" Done.")

    # --- Prepare to output video (?) -----------------------------------------
    if opts.capture_directory != None:
        from signal import SIG_DFL, signal

        # SIGPIPE is not used on any Windows system.
        if not sys.platform.startswith('win'):
            from signal import SIGPIPE
        base_dir = opts.capture_directory
        MOVIE_DIRECTORY = base_dir
        DEBUG_DIRECTORY = os.path.join(base_dir, DEBUG_SUBDIRECTORY)
        FRAME_DIRECTORY = os.path.join(base_dir, FRAME_SUBDIRECTORY)
        for d in [base_dir, MOVIE_DIRECTORY, DEBUG_DIRECTORY, FRAME_DIRECTORY]:
            if not os.path.isdir(d):
                os.mkdir(d)
        os.environ["SDL_VIDEODRIVER"] = "dummy"
        if opts.debug:
            # TODO: maybe opts.verbose
            print("SDL_VIDEODRIVER set to 'dummy'")
    else:
        FRAME_DIRECTORY = ""

    # --- Initialize Surface State --------------------------------------------
    if init_state:
        grid = init_state
    else:
        if opts.grid is None:
            raise Exception("Initial grid state required.")
        grid = opts.grid

    # --- Setup for either Asynchronous or Synchronous Simulation -------------
    if opts.simulation_type == "asynchronous":
        if opts.debug:
            print("Grid is type " + str(type(grid)))
            print("Initializing simulator with surface:\n" + str(grid))
            for x in range(grid.x_size):
                for y in range(grid.y_size):
                    print("(" + str(x) + "," + str(y) + "): " + str(grid.grid[x, y]))
        simulation = QueueSimulator(surface=grid,
                                    transition_rules=opts.transition_rules,
                                    seed=opts.rng_seed,
                                    group_selection_seed=group_selection_seed,
                                    simulation_duration=opts.max_duration,
                                    rxns=rxns
                                    )
        simulation.init_wall_time = process_time()
    # TODO: support synchronous mode.
    elif opts.simulation_type == "synchronous":
        # TODO: study when to use it;
        # TODO: update the synchronous simulator to be the same fashion as the queue simulator
        simulation = SynchronousSimulator(
            surface=grid,
            update_rule=opts.update_rule,
            seed=opts.rng_seed,
            simulation_duration=opts.max_duration,
        )
        simulation.init_wall_time = process_time()
    else:
        raise Exception('Unknown simulation type "' + opts.simulation_type + '".')
    time = simulation.time
    event_history = EventHistory()
    # TODO(Andrew): These next four lines are for integration with our stuff
    simulation.rxns = rxns

    # Properties for visualizing the spectra along the grid frames.
    simulation.running_average = running_average
    simulation.spectra_in_video = spectra_in_video
    simulation.ir_intensities = ir_intensities

    ################
    # PYGAME SETUP #
    ################
    if opts.debug:
        print("Beginning Pygame setup...")

    # --- Get Display Class from Default Grid Types ---------------------------
    if opts.grid_type == 'parallel_emulated':
        grid_display = ParallelEmulatedSquareGridDisplay(grid=grid,
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
                                                         display_text=opts.display_text)
    elif opts.grid_type == 'standard':
        if display_class:
            DisplayClass = display_class
        elif opts.grid_type == "standard" and opts.surface_geometry == "square":
            DisplayClass = SquareGridDisplay
        elif opts.grid_type == "standard" and opts.surface_geometry == "hex":
            DisplayClass = HexGridDisplay
        grid_display = DisplayClass(grid=grid,
                                    colormap=opts.COLORMAP,
                                    min_x=MIN_GRID_WIDTH,
                                    min_y=0,
                                    pixels_per_node=opts.pixels_per_node,
                                    display_text=opts.display_text)
    else:
        raise Exception("Unrecognized grid type '" + opts.grid_type + "'")

    # --- Create Displays for Legend, Title, etc. -----------------------------
    legend_display = LegendDisplay(colormap=opts.COLORMAP)

    # Width only requires legend and grid sizes to calculate
    display_width = grid_display.display_width + legend_display.display_width

    # Width used to calculate time label and button placements
    time_display = TimeDisplay(display_width)

    # Display for the additional title
    title_display = TextDisplay(display_width, text="Surface CRN Trajectory")
    show_title = bool(spectra_in_video)

    # --- Create Buttons for Synchronous Simulation ---------------------------
    button_y = time_display.display_height + grid_display.display_height + 1
    # max(legend_display.display_height, grid_display.display_height) + 1
    # (int(display_width/2) - (button_width + button_buffer), button_y,
    play_back_button = PygButton(rect=
                                 (legend_display.display_width + button_buffer, button_y,
                                  button_width, button_height),
                                 caption='<<')
    step_back_button = PygButton(rect=
                                 (play_back_button.rect.right + button_buffer, button_y,
                                  button_width, button_height),
                                 caption='< (1)')
    pause_button = PygButton(rect=
                             (step_back_button.rect.right + button_buffer, button_y,
                              button_width, button_height),
                             caption='Pause')
    step_button = PygButton(rect=
                            (pause_button.rect.right + button_buffer, button_y,
                             button_width, button_height),
                            caption='(1) >')
    play_button = PygButton(rect=
                            (step_button.rect.right + button_buffer, button_y,
                             button_width, button_height),
                            caption='>>')
    clip_button = PygButton(rect=
                            (play_button.rect.right + 4 * button_buffer, button_y,
                             button_width * 1.1, button_height),
                            caption='Uncache')

    if show_title:
        display_height = max(legend_display.display_height + 2 * legend_display.VERTICAL_BUFFER +
                             time_display.display_height + title_display.display_height,
                             button_y + button_height + 2 * button_buffer)
    else:
        display_height = max(legend_display.display_height + \
                             2 * legend_display.VERTICAL_BUFFER + time_display.display_height,
                             button_y + button_height + 2 * button_buffer)

    # --- Initialize PyGame Canvas --------------------------------------------
    if opts.debug:
        print("Initializing display of size " + str(display_width) + ", " +
              str(display_height) + ".")
    display_surface = pygame.display.set_mode((display_width,
                                               display_height), 0, 32)
    simulation.display_surface = display_surface
    # except:
    #     display_surface = pygame.display.set_mode((display_width,
    #                                                display_height))
    if opts.debug:
        print("Display initialized. Setting caption.")
    pygame.display.set_caption('Surface CRN Simulator')
    if opts.debug:
        print("Caption set, initializing clock.")
    fpsClock = pygame.time.Clock()

    if opts.debug:
        print("Clock initialized. Filling display with white.")
    # Initial render
    display_surface.fill(WHITE)
    if opts.debug:
        # TODO: maybe this should be shown in verbose mode
        print("Pygame setup done, first render attempted.")

    # Make the options menu.
    # opts_menu = MainOptionMenu()
    # opts_menu.update()

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
    #  TODO: - 200 is a trial for x_pos

    # legend_display.display_height is intended to be a very small number.
    grid_display.render(display_surface, x_pos=legend_display.display_width,
                        y_pos=time_display.y_pos + time_display.display_height,
                        width=time_display.display_width - legend_display.display_width,
                        height=legend_display.display_height * 4)
    # grid_display.render(display_surface, x_pos=time_display.x_pos,
    #                                      y_pos=time_display.y_pos +
    #                                              time_display.display_height
    #                                      )

    if opts.saving_movie:
        simulation.display_surface_size = display_height * display_width
    else:
        play_back_button.draw(display_surface)
        step_back_button.draw(display_surface)
        pause_button.draw(display_surface)
        step_button.draw(display_surface)
        play_button.draw(display_surface)
        clip_button.draw(display_surface)

    pygame.display.flip()
    # TODO
    # progress_bar = ProgressBar(total_tasks=round(opts.max_duration/opts.fps))
    progress_bar = None
    update_display(opts, simulation, progress_bar, grid_display, FRAME_DIRECTORY,
                   time_display=time_display,
                   title_display=title_display,
                   spectra_max_conc=spectra_max_conc)

    # State variables for simulation
    next_reaction_time = 0
    prev_reaction_time = 0
    next_reaction = None
    prev_reaction = None
    running = True  # TEMPORARY FIX ME!!!
    first_frame = True
    last_frame = False
    running_backward = False

    # --- Run PyGame Loop, i.e. actually Simulate -----------------------------
    if opts.debug:
        print("Beginning simulation....")
    # Iterate through events
    while True:
        # Check for interface events
        for event in pygame.event.get():
            if 'click' in play_back_button.handleEvent(event):
                running = True
                running_backward = True
                last_frame = False
            if 'click' in step_back_button.handleEvent(event):
                running = False
                last_frame = False
                if event_history.at_beginning():
                    time = 0
                    time_display.time = 0
                    time_display.render(display_surface, x_pos=time_display.x_pos, y_pos=time_display.y_pos)
                    pygame.display.update()
                else:
                    prev_reaction = event_history.previous_event()
                    event_history.increment_event(-1)
                    prev_reaction_rule = prev_reaction.rule
                    for i in range(len(prev_reaction.participants)):
                        cell = prev_reaction.participants[i]
                        state = prev_reaction_rule.inputs[i]
                        cell.state = state
                    display_next_event(prev_reaction, grid_display)
                    if event_history.at_beginning():
                        time = 0
                    else:
                        # Note that this is NOT the same as the time from the
                        # "previous event" that we're currently processing --
                        # it's actually the time of the event *before* the one
                        # we just undid.
                        time = event_history.previous_event().time
                    time_display.time = time
                    time_display.render(display_surface, x_pos=time_display.x_pos, y_pos=time_display.y_pos)
                    pygame.display.update()
            if 'click' in pause_button.handleEvent(event):
                running = False
            if 'click' in step_button.handleEvent(event):
                running = False
                first_frame = False
                # Process a single reaction
                reading_history = False
                if not event_history.at_end():
                    reading_history = True
                    next_reaction = event_history.next_event()
                    event_history.increment_event(1)
                    next_reaction_rule = next_reaction.rule
                    for i in range(len(next_reaction.participants)):
                        cell = next_reaction.participants[i]
                        state = next_reaction_rule.outputs[i]
                        cell.state = state
                if not next_reaction:
                    next_reaction = simulation.process_next_reaction()
                    if not reading_history:
                        event_history.add_event(next_reaction)
                        event_history.increment_event(1)
                next_reaction_time = next_reaction.time
                display_next_event(next_reaction, grid_display)
                time = next_reaction_time
                time_display.time = time
                time_display.render(display_surface, x_pos=time_display.x_pos, y_pos=time_display.y_pos)
                pygame.display.update()
                next_reaction = None
                if opts.debug:
                    print("State after update: " + str(grid))
            if 'click' in play_button.handleEvent(event):
                running = True
                running_backward = False
                first_frame = False
            if 'click' in clip_button.handleEvent(event):
                event_history.clip()
                simulation.time = time
                simulation.reset()
                if opts.debug:
                    for rxn in list(simulation.event_queue.queue):
                        print(rxn)
            if event.type == QUIT:
                if opts.saving_movie:
                    movie_file.close()
                cleanup_and_exit(simulation)
        # Don't do anything if paused.
        if not running:
            pygame.display.update()
            continue

        # Update time
        if opts.debug:
            print(f"Updating time: time = {time}, running_backward = "
                  f"{running_backward}, first_frame = {first_frame}, "
                  f"last_frame = {last_frame}")
        if running_backward and not first_frame:
            # prev_reaction_time = time
            time -= opts.speedup_factor * 1. / opts.fps
            last_frame = False
        elif not running_backward and not last_frame:
            # next_reaction_time = time
            time += opts.speedup_factor * 1. / opts.fps
            first_frame = False
        if opts.debug:
            print(f"Updating time to {time}")
        time_display.time = time
        time_display.render(display_surface, x_pos=time_display.x_pos, y_pos=time_display.y_pos)

        # Process any simulation events that have happened since the last tick.
        if opts.debug:
            print("Checking for new events...")
        if running_backward and not first_frame:
            if opts.debug and not event_history.at_beginning():
                print(f"While running backwards, checking if there are any "
                      f"events: time = {time}, previous event time = "
                      f"{event_history.previous_event().time}")
            while not event_history.at_beginning() and \
                    event_history.previous_event().time > time:
                prev_reaction = event_history.previous_event()
                previous_reaction_time = prev_reaction.time
                if event_history.at_beginning():
                    first_frame = True
                event_history.increment_event(-1)
                # if opts.debug:
                #     print("While running backwards, undoing reaction: "
                #           f"{prev_reaction}")
                for i in range(len(prev_reaction.participants)):
                    cell = prev_reaction.participants[i]
                    state = prev_reaction.rule.inputs[i]
                    cell.state = state
                next_reaction_time = prev_reaction_time
                prev_reaction_time = prev_reaction.time if prev_reaction \
                    else 0
                if opts.debug:
                    print("Displaying a new event")
                display_next_event(prev_reaction, grid_display)
                if opts.debug and not event_history.at_beginning():
                    print(f"While running backwards, checking if there are any "
                          f"events: time = {time}, previous event time = "
                          f"{event_history.previous_event().time}")
        elif not running_backward and not last_frame:
            while (not event_history.at_end() or not simulation.done()) \
                    and next_reaction_time < time:
                if event_history.at_end():
                    # Advance the reaction
                    next_reaction = simulation.process_next_reaction()
                    if next_reaction:
                        event_history.add_event(next_reaction)
                        event_history.increment_event(1)
                else:
                    print("reading from event history")
                    next_reaction = event_history.next_event()
                    event_history.increment_event(1)
                    # TODO: does this update the visuals?
                    for i in range(len(next_reaction.participants)):
                        cell = next_reaction.participants[i]
                        state = next_reaction.rule.outputs[i]
                        cell.state = state

                prev_reaction_time = next_reaction_time
                next_reaction_time = next_reaction.time if next_reaction \
                    else opts.max_duration + 1
                if opts.debug:
                    print("Displaying a new event")
                display_next_event(next_reaction, grid_display)

        # Render updates and make the next clock tick.
        if opts.debug:
            print("Updating display.")
        update_display(opts, simulation, progress_bar, grid_display, FRAME_DIRECTORY, time_display=time_display,
                       title_display=title_display,
                       spectra_max_conc=spectra_max_conc)
        fpsClock.tick(opts.fps)

        # Check for simulation completion...
        if opts.debug:
            print("Checking for simulation completion...")
        if event_history.at_end() and running and not running_backward and \
                (simulation.done() or time > opts.max_duration):
            if opts.debug:
                print("Done! Cleaning up now.")
            last_frame = True
            running = False
            # Set the time to final time when done.
            time = time
            # time = opts.max_duration
            time_display.time = time
            time_display.render(display_surface, x_pos=time_display.x_pos,
                                y_pos=time_display.y_pos)  # opts_menu.display_height)
            if next_reaction:
                display_next_event(next_reaction, grid_display)
            update_display(opts, simulation, progress_bar, grid_display, FRAME_DIRECTORY, time_display=time_display,
                           title_display=title_display,
                           spectra_max_conc=spectra_max_conc, check_terminate=True)
            if opts.debug:
                print("Simulation state at final time " + \
                      str(opts.max_duration) + ":")
                print(str(grid))
            if opts.capture_directory != None:
                # Use ffmpeg to convert images to movie.
                if os.name == 'nt':
                    ffmpeg_name = 'ffmpeg.exe'
                elif os.name == 'posix':
                    # ffmpeg_name = 'ffmpeg'
                    name_found = False
                    for possible_name in ['/usr/local/bin/ffmpeg',
                                          '/usr/bin/ffmpeg']:
                        if os.path.isfile(possible_name):
                            ffmpeg_name = possible_name  # EW: not finding it?
                            name_found = True
                            break
                    if not name_found:
                        ffmpeg_name = 'ffmpeg'
                        print("Could not find executable ffmpeg in" +
                              " any of the expected locations!")
                        # raise Exception("Could not find executable ffmpeg in"
                        #                 " any of the expected locations!")
                else:
                    raise Exception("Unexpected OS name '" + os.name + "'")

                if not sys.platform.startswith('win'):
                    signal(SIGPIPE, SIG_DFL)
                width = display_surface.get_width()
                height = display_surface.get_height()
                movie_filename = os.path.join("", MOVIE_DIRECTORY,
                                              opts.movie_title + ".mp4")
                if opts.debug:
                    print("Writing movie to  file " + movie_filename +
                          "\n")
                command = [ffmpeg_name,
                           '-y',  # Overwrite output file
                           # '-framerate', str(opts.fps),
                           '-start_number', '1',  # EW: might it default to
                           # starting with 0?
                           '-i', os.path.join(FRAME_DIRECTORY,
                                              opts.movie_title + "_%d.jpeg"),
                           '-an',  # no audio
                           # Width and height need to be divisible by 2.
                           # Round up if necessary.
                           '-vf', 'pad=ceil(iw/2)*2:ceil(ih/2)*2',
                           movie_filename
                           ]

                if opts.debug:
                    print("Calling ffmpeg with: " + str(command))
                    print("And right now the current dir is " + os.getcwd())
                    print("opts.capture_directory = " + opts.capture_directory)
                    # if opts.debug:

                    print("Writing movie with command:\n")
                    print("\t" + str(command) + "\n")
                debug_output_stream = open(os.path.join(opts.capture_directory,
                                                        "debug",
                                                        "ffmpeg_debug.dbg"), 'w')
                print('command', command, '\n')
                print('dos', debug_output_stream, '\n')
                proc = sp.Popen(command,
                                stdout=debug_output_stream,
                                stderr=sp.STDOUT)
                proc.communicate()
                if opts.debug:
                    print("Finished ffmpeg call.")

                return simulation.concs, simulation.times
        if event_history.at_beginning() or time == 0:
            first_frame = True


# end def main()

def display_next_event(next_reaction, grid_display):
    # Update the display for a reaction. Reaction will be an Event object (see
    # surface_crns.simulators.py).
    DEBUG = False

    if DEBUG:
        print("Moving to next reaction:")

    if not next_reaction:
        if DEBUG:
            print("No reaction returned from event queue. Finished")
        return -1

    next_reaction_time = next_reaction.time

    if DEBUG:
        print("Updating display based on event " +
              str(next_reaction))

    # Display any changes made
    participants = next_reaction.participants
    inputs = next_reaction.rule.inputs
    outputs = next_reaction.rule.outputs
    # Update reactants (if changed)
    for i in range(len(participants)):
        # TODO: this is currently very expensive for coord_grid's Voronoi pictures.
        # TODO: accomodations for 2-sized species
        if i > len(inputs) - 1:
            grid_display.update_node(participants[i])
        elif inputs[i] != outputs[i]:
            grid_display.update_node(participants[i])
        elif DEBUG:
            print("Input " + str(i + 1) + " and output " + str(i + 1) + " match " +
                  "for rule " + str(next_reaction.rule) + "\ncell " +
                  str(participants[0].position) + " not updated.")

    # Display current state to stdout.
    if DEBUG:
        print("Simulation state at time " + str(next_reaction_time) + \
              ":\n" + str(grid_display.grid))

    return next_reaction_time


def cleanup_and_exit(simulation):
    pygame.quit()
    print("Simulation state at termination (T = " + str(simulation.time) + "):")
    print(str(simulation.surface))
    sys.exit()


def update_display(opts, simulation, progress_bar, grid_display, FRAME_DIRECTORY=None, time_display=None,
                   title_display=None,
                   spectra_max_conc=-1, check_terminate=False):
    if simulation.times:
        time = simulation.times[-1]
    else:
        time = 0
    if check_terminate:
        if math.isclose(time, opts.max_duration, abs_tol=0.1):
            text = "Generating video frames completed"
        else:
            # TODO: don't let this take 2 lines.
            text = f"Generating video frames completed; simulation terminated early at {time:.1f} s"

        ipython_visuals.update_progress(time / opts.max_duration, text,
                                        beginning=time == 0,
                                        terminating=True)
    else:

        ipython_visuals.update_progress(time / opts.max_duration, "Generating video frames",
                                        beginning=time == 0,
                                        terminating=check_terminate)
    # TODO: display a message on termination
    # and
    # progress_bar.bar()

    grid_display.re_render()
    # TODO
    # print(type(simulation))
    if opts.capture_directory == None:
        pygame.display.update()
        pygame.display.flip()
    else:
        if opts.debug:
            print("capture directory is: " + str(opts.capture_directory))
        if FRAME_DIRECTORY == None:
            raise Exception("FRAME_DIRECTORY should be set if a capture" +
                            " directory is set.")
        try:
            capture_time = simulation.capture_time
        except AttributeError:
            simulation.capture_time = 0
            capture_time = 0

        try:
            frame_number = simulation.frame_number
        except AttributeError:
            simulation.frame_number = 1
            frame_number = 1

        if simulation.time >= capture_time:
            if opts.debug:
                print("movie title is: " + str(opts.movie_title))
            frame_filename = os.path.join(FRAME_DIRECTORY, opts.movie_title
                                          + "_" + str(frame_number) +
                                          ".jpeg")
            if opts.debug:
                print("Saving frame at: " + frame_filename)

            # TODO: create and save more contents:
            # Currently this block adds the Gaussian figures
            screen = simulation.display_surface
            if simulation.spectra_in_video:
                trajectory_size = simulation.display_surface.get_size()
                trajectory_width = trajectory_size[0]
                trajectory_height = trajectory_size[1]
                # Horizontal gap between grid video and spectrum.
                h_gap = 40
                up_gap = 50

                # Make entire screen, display the saved video on the left side.
                spectrum_width = trajectory_width * 1 / 2
                # The number of spectra in display.
                number_of_spectra = 3
                # The number of spectra per column
                number_of_spectrum_rows = 2
                # The number of columns of spectra
                number_of_spectrum_cols = number_of_spectra // number_of_spectrum_rows \
                    if number_of_spectra % number_of_spectrum_rows == 0 else \
                    number_of_spectra // number_of_spectrum_rows + 1

                entire_screen = pygame.Surface([trajectory_width + (spectrum_width + h_gap) * number_of_spectrum_cols,
                                                trajectory_height])
                entire_screen.fill((255, 255, 255))
                entire_screen.blit(simulation.display_surface, (0, 0))

                r = Results.from_counts(simulation.rxns, simulation.surface.species_count())
                dpi = 100

                if time_display is not None:
                    title_x, title_y = title_display.x_pos, title_display.y_pos

                    display = TextDisplay(spectrum_width, text="Dynamic XPS Spectrum")
                    display.render(entire_screen, title_x + trajectory_width + h_gap, title_y)

                    # TODO: add a text display

                    if time_display.get_time() > simulation.running_average:
                        start_time = time_display.get_time() - simulation.running_average
                    else:
                        start_time = 0
                    if time_display.get_time() == 0:
                        time_period_string = f"T = {time_display.get_time():.2f}"
                    else:
                        time_period_string = f"T = {start_time:.2f} to T = {time_display.get_time():.2f}"
                    running_avg_display = TextDisplay(spectrum_width, font_size=18,
                                                      text=time_period_string)

                    # Don't save time display's x, y locations.
                    time_x, time_y = time_display.x_pos, time_display.y_pos
                    time = time_display.get_time()
                    new_time_display = TimeDisplay(spectrum_width)
                    new_time_display.set_time(time)
                    new_time_display.render(entire_screen, x_pos=time_x + trajectory_width + h_gap, y_pos=time_y)
                    gap = 0
                    fig_gap = 0  # the gap between two figures
                    fig_height = trajectory_height - title_display.display_height - \
                                 new_time_display.display_height - running_avg_display.display_height - gap - fig_gap

                    # Fit two pictures
                    fig_height = fig_height / 2
                    y_lim = round(1.1 * spectra_max_conc) if spectra_max_conc != -1 else simulation.surface.num_nodes
                    # TODO: also set an X axis limit
                    raw_data, size = r.raw_string_gaussian(y_upper_limit=y_lim,
                                                           fig_size=(spectrum_width / dpi, fig_height / dpi),
                                                           dpi=dpi)
                    gaussian = pygame.image.fromstring(raw_data, size, "RGB")

                    start = new_time_display.y_pos + new_time_display.display_height + gap

                    # Blit the entire Gaussian into certain positions on the screen.
                    entire_screen.blit(gaussian, (new_time_display.x_pos, start))

                    start += fig_height
                    running_avg_display.render(entire_screen, x_pos=new_time_display.x_pos, y_pos=start)

                    if simulation.concs:
                        r = Results.from_concs_times(None, simulation.rxns, simulation.concs, simulation.times)
                    else:
                        r = r

                    # print(r.df_raw)
                    if simulation.time > r.df_raw.index.max():
                        simulation_time = r.df_raw.index.max()
                    else:
                        simulation_time = simulation.time

                    starting_time = max(0, simulation_time - simulation.running_average)
                    # print("calculating running average")
                    # print("starting time", starting_time)
                    # print("duration", simulation.running_average)

                    # TODO: neatly deal with early termination.
                    # print("sim time", simulation_time)

                    raw_data, size = r.raw_string_gaussian(y_upper_limit=y_lim,
                                                           t=starting_time,
                                                           avg_duration=simulation.running_average,
                                                           fig_size=(spectrum_width / dpi, fig_height / dpi),
                                                           dpi=dpi)
                    gaussian = pygame.image.fromstring(raw_data, size, "RGB")
                    start += running_avg_display.display_height + fig_gap
                    entire_screen.blit(gaussian, (new_time_display.x_pos, start))

                    # Display IR
                    # Cascading X positions from the previous display
                    left_x_pos = new_time_display.x_pos + spectrum_width + h_gap

                    display = TextDisplay(spectrum_width, text="Dynamic Concentration Plot")
                    display.render(entire_screen, left_x_pos, title_y)
                    time = new_time_display.get_time()
                    time_display = TimeDisplay(spectrum_width)
                    time_display.set_time(time)
                    time_display.render(entire_screen, x_pos=left_x_pos, y_pos=new_time_display.y_pos)

                    # TODO(Andrew): This is what generates that warning about
                    #  xlim and left == right == 0 and singular tranformation.
                    raw_string_c, size = r.raw_string_evolution(ylim=y_lim,
                                                                return_fig=True,
                                                                fig_size=(spectrum_width / dpi, fig_height / dpi),
                                                                dpi=dpi,)

                    gaussian = pygame.image.fromstring(raw_string_c, size, "RGB")

                    entire_screen.blit(gaussian, (left_x_pos,
                                                  new_time_display.y_pos + new_time_display.display_height + gap))

                    # Plot the rolling mean below the normal one
                    def rolling_mean(df, window):
                        def rolled_row(row):
                            return df[max(0, time - window):row.name].mean()
                        return rolled_row
                    # r.df_raw.apply(rolling_mean(r.df_raw, 1), axis=1)

                    raw_string_c, size = r.raw_string_evolution(ylim=y_lim,
                                                                return_fig=True,
                                                                fig_size=(spectrum_width / dpi, fig_height / dpi),
                                                                dpi=dpi,
                                                                use_raw_data=False)

                    evolution = pygame.image.fromstring(raw_string_c, size, "RGB")
                    entire_screen.blit(evolution, (left_x_pos, start))

                    running_avg_display = TextDisplay(spectrum_width, font_size=18,
                                                      text=time_period_string)
                    running_avg_display.render(entire_screen,
                                               x_pos=left_x_pos,
                                               y_pos=(
                            start - new_time_display.display_height - gap))

                else:
                    raw_data, size = r.raw_string_gaussian(y_upper_limit=simulation.surface.num_nodes,
                                                           fig_size=(trajectory_width / dpi,
                                                                     (trajectory_height - up_gap) / dpi),
                                                           dpi=dpi)
                    gaussian = pygame.image.fromstring(raw_data, size, "RGB")
                    entire_screen.blit(gaussian, (trajectory_width + h_gap, up_gap - 10))

                screen = entire_screen
            pygame.image.save(screen, frame_filename)

            # Determine next capture time
            simulation.capture_time = capture_time + 1. / opts.capture_rate
            simulation.frame_number = frame_number + 1

            # Check the space used. If it's too much, save one last frame
            # and terminate.
            try:
                simulation.pixels_saved += simulation.display_surface_size
            except AttributeError:
                simulation.pixels_saved = simulation.display_surface_size

            terminate = False
            if simulation.pixels_saved > CUTOFF_SIZE:
                termination_string = "Simulation terminated at maximum amount of video frame pixels saved " + \
                                     str(simulation.pixels_saved) + \
                                     " pixels (~" + \
                                     str(CUTOFF_SIZE / 10000000) + " Mb)."
                teriminate = True

            # Check the timer. If it's been more than an hour, terminate.
            if process_time() - simulation.init_wall_time > CUTOFF_TIME:
                termination_string = "Simulation terminated at maximum allowed " \
                                     f"processing time t = {CUTOFF_TIME}s"
                terminate = True

            if terminate:
                # TODO: clean up this section
                print(termination_string)

                display_width = simulation.display_surface.get_size()[1]
                text_display = TextDisplay(display_width)
                text_display.text = termination_string
                text_display.render(simulation.display_surface, x_pos=0, y_pos=0)
                frame_filename = os.path.join(FRAME_DIRECTORY,
                                              opts.movie_title + "_" +
                                              str(frame_number) + ".jpeg")
                if opts.debug:
                    print("Saving final frame at: " + frame_filename)
                pygame.image.save(simulation.display_surface, frame_filename)

                cleanup_and_exit(simulation)


if __name__ == '__main__':
    if PROFILE:
        """try:
            import statprof
            statprof.start()
            try:
                main()
            finally:
                statprof.stop()
                statprof.display()
        except ImportError:"""
        cProfile.run("main()", sort='tottime')
    else:
        main()
