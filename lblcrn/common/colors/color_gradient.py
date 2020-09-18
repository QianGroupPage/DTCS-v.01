# based on https://bsou.io/posts/color-gradients-with-python

from matplotlib import colors

def hex_to_RGB(hex):
  ''' "#FFFFFF" -> [255,255,255] '''
  # Pass 16 to the integer function for change of base
  return [int(hex[i:i+2], 16) for i in range(1,6,2)]


def RGB_to_hex(RGB):
  ''' [255,255,255] -> "#FFFFFF" '''
  # Components need to be integers for hex to make sense
  RGB = [int(x) for x in RGB]
  return "#"+"".join(["0{0:x}".format(v) if v < 16 else
            "{0:x}".format(v) for v in RGB])


def color_to_RGB(color):
    # TODO: round first 3 values if list or tuple
    if isinstance(color, (tuple, list)):
        color = tuple(round(color) for color in color[:3])
        return color
    color = colors.to_rgb(color)
    color = tuple(round(val * 255) for val in color)
    return color


def color_to_HEX(color, zero_to_one_range=True):
    # Send to matplotlib.color early if termination is early.
    for c in color:
        if isinstance(c, str):
            return tuple(colors.to_hex(color))
    if not zero_to_one_range:
        color = tuple(c / 256 for c in color)
    return tuple(colors.to_hex(color))


def color_dict(gradient):
  ''' Takes in a list of RGB sub-lists and returns dictionary of
    colors in RGB and hex form for use in a graphing function
    defined later on '''
  return {"hex":[RGB_to_hex(RGB) for RGB in gradient],
      "r":[RGB[0] for RGB in gradient],
      "g":[RGB[1] for RGB in gradient],
      "b":[RGB[2] for RGB in gradient]}


def linear_gradient(start_hex, finish_hex="#FFFFFF", n=10):
  ''' returns a gradient list of (n) colors between
    two hex colors. start_hex and finish_hex
    should be the full six-digit color string,
    inlcuding the number sign ("#FFFFFF") '''
  # Starting and ending colors in RGB form
  s = hex_to_RGB(start_hex)
  f = hex_to_RGB(finish_hex)
  # Initilize a list of the output colors with the starting color
  RGB_list = [s]
  # Calcuate a color at each evenly spaced value of t from 1 to n
  for t in range(1, n):
    # Interpolate RGB vector for color at the current value of t
    curr_vector = [
      int(s[j] + (float(t)/(n-1))*(f[j]-s[j]))
      for j in range(3)
    ]
    # Add it to our list of output colors
    RGB_list.append(curr_vector)

  return color_dict(RGB_list)

def color_series(start_hex, n=10, finish_hex="#FFFFFF"):
    # Lighter colors are shown first
    return linear_gradient(start_hex, finish_hex=finish_hex, n=n + 1)["hex"][:-1][::-1]