from IPython.display import clear_output


def update_progress(progress, header="Progress"):
    # https://www.mikulskibartosz.name/how-to-display-a-progress-bar-in-jupyter-notebook/
    bar_length = 300
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



