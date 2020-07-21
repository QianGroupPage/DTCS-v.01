from IPython.display import clear_output
import os

def update_progress(progress, header="Progress", beginning=False, terminating=False):
    # https://www.mikulskibartosz.name/how-to-display-a-progress-bar-in-jupyter-notebook/
    bar_length = 80
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
    if progress < 0:
        progress = 0
    if progress >= 1:
        progress = 1

    block = int(round(bar_length * progress))
    clear_output(wait=True)

    pound_dashes = "#" * block + "-" * (bar_length - block)
    text = f"{header}: [{pound_dashes}] {progress * 100:.1f}%"
    print(text)

    terminal_output = open('/dev/stdout', 'w')
    terminal_columns = os.get_terminal_size().columns
    if len(text) > terminal_columns:
        bar_length = bar_length - (len(text) - terminal_columns)
        pound_dashes = "#" * block + "-" * (bar_length - block)
        text = f"{header}: [{pound_dashes}] {progress * 100:.1f}%"
    if not terminating:
        terminal_output.write(f"{text} \r")
        terminal_output.flush()
    else:
        terminal_output.write(f"{text}\n")

    # Prevent the progress bar from showing on command line
    # sys.stdout.flush()


class ProgressBar:
    def __init__(self, header="Progress", collect_time_stamps=False, total_tasks=1):
        pass
    #     self.header = header
    #     self.collect_time_stamps = False
    #     if self.collect_time_stamps:
    #         self.time_stamps = []
    #     self.bar = alive_bar(total_tasks)
    #
    # def progress_and_print(self, progress):
    #     update_progress(progress, self.header)
    #
    # def bar(self):
    #     clear_output(wait=True)
    #     self.bar()



