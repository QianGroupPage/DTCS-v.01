
from IPython.core import interactiveshell


def format_workaround(obj, include=None, exclude=None):
    """This is word-for-word copied from InteractiveShell.format, but without
    the part that checks for _ipython_display_. This allows our Monty objects
    to still display even when Crystal Toolkit adds an _ipython_display_ to
    them. If they have nothing else to show, they'll use the Crystal Toolkit
    version. This should eventually be removed; I think Crystal Toolkit should
    add a type of formatter, or something. Maybe I'll fork it and give it a
    go someday."""
    fmt = interactiveshell.InteractiveShell.instance().display_formatter

    format_dict, md_dict = fmt.mimebundle_formatter(obj, include=include, exclude=exclude)

    if format_dict or md_dict:
        if include:
            format_dict = {k: v for k, v in format_dict.items() if k in include}
            md_dict = {k: v for k, v in md_dict.items() if k in include}
        if exclude:
            format_dict = {k: v for k, v in format_dict.items() if k not in exclude}
            md_dict = {k: v for k, v in md_dict.items() if k not in exclude}

    for format_type, formatter in fmt.formatters.items():
        if format_type in format_dict:
            # already got it from mimebundle, maybe don't render again.
            # exception: manually registered per-mime renderer
            # check priority:
            # 1. user-registered per-mime formatter
            # 2. mime-bundle (user-registered or repr method)
            # 3. default per-mime formatter (e.g. repr method)
            try:
                formatter.lookup(obj)
            except KeyError:
                # no special formatter, use mime-bundle-provided value
                continue
        if include and format_type not in include:
            continue
        if exclude and format_type in exclude:
            continue

        md = None
        try:
            data = formatter(obj)
        except:
            # FIXME: log the exception
            raise

        # formatters can return raw data or (data, metadata)
        if isinstance(data, tuple) and len(data) == 2:
            data, md = data

        if data is not None:
            format_dict[format_type] = data
        if md is not None:
            md_dict[format_type] = md

    return format_dict, md_dict
