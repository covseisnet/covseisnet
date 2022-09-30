# coding: utf-8

import os

import numpy as np
import matplotlib.pyplot as plt

from scipy import signal

SIZE_FIGURE_INCHES = 2, 2
SIZE_FIGURE_SMALL_INCHES = 1, 1
SIZE_TRIANGLES = 15
SIZE_TEXT = 22

WAVE_LINEWIDTH = 3
WAVE_STRETCH_X = 0.8
WAVE_STRETCH_Y = 0.7
WAVE_FREQUENCY = 5
WAVE_WIDTH = 6

BRACKETS_LINEWIDTH = 4
BRACKET_HORIZONTAL_MARGIN = 0.75
BRACKET_VERTICAL_MARGIN = 0.75

SAVEPATH = os.path.dirname(__file__)
FILEPATH_LOGO_1 = os.path.join(SAVEPATH, "logo_square.svg")
FILEPATH_LOGO_2 = os.path.join(SAVEPATH, "logo_square_small.svg")

COLOR_1 = "#4B8BBE"
COLOR_2 = "#FFD43B"
COLOR_3 = "0.5"

plt.rcParams["figure.dpi"] = 200
plt.rcParams["savefig.transparent"] = True
plt.rcParams["savefig.bbox"] = "tight"
plt.rcParams["lines.dash_joinstyle"] = "round"


def item_coordinates():
    coordinates = np.meshgrid(np.arange(3), np.arange(3))
    coordinates = coordinates[0].ravel()[1:], coordinates[1].ravel()[1:]
    return coordinates


def item_coordinates_small():
    coordinates = np.meshgrid(np.arange(2) * 2, np.arange(2) * 2)
    coordinates = coordinates[0].ravel()[1:], coordinates[1].ravel()[1:]
    return coordinates


def logo_square():

    # Figure and axes instance
    figure, axes = plt.subplots(1, figsize=SIZE_FIGURE_INCHES)

    # Plot triangles
    for coordinates in zip(*item_coordinates()):
        x, y = coordinates
        if 2 - x == y:
            color = COLOR_1
        elif (x == 2) & (y == 2):
            color = COLOR_3
        else:
            color = COLOR_2
        axes.plot(
            *coordinates, "^", ms=SIZE_TRIANGLES, clip_on=False, color=color
        )

    # Plot square braquets
    hspace = BRACKET_HORIZONTAL_MARGIN
    wspace = BRACKET_VERTICAL_MARGIN

    # Bottom left
    plt.plot(
        [-hspace, 0],
        [-wspace, -wspace],
        color=COLOR_3,
        linewidth=BRACKETS_LINEWIDTH,
    )

    # Bottom right
    plt.plot(
        [2 + hspace, 2],
        [-wspace, -wspace],
        color=COLOR_3,
        linewidth=BRACKETS_LINEWIDTH,
    )

    # Top left
    plt.plot(
        [-hspace, 0],
        [2 + wspace, 2 + wspace],
        color=COLOR_3,
        linewidth=BRACKETS_LINEWIDTH,
    )

    # Top right
    plt.plot(
        [2, 2 + hspace],
        [2 + wspace, 2 + wspace],
        color=COLOR_3,
        linewidth=BRACKETS_LINEWIDTH,
    )

    # Left
    plt.plot(
        [-hspace, -hspace],
        [-wspace, 2 + wspace],
        color=COLOR_3,
        linewidth=BRACKETS_LINEWIDTH,
    )

    # Right
    plt.plot(
        [2 + hspace, 2 + hspace],
        [-wspace, 2 + wspace],
        color=COLOR_3,
        linewidth=BRACKETS_LINEWIDTH,
    )

    # Text
    plt.text(
        1,
        -1,
        "covseisnet",
        name="Ubuntu",
        va="top",
        ha="center",
        size=SIZE_TEXT,
        color=COLOR_3,
        weight=800,
    )
    # Wave
    time = np.linspace(-0.5, 0.5) * WAVE_STRETCH_X
    wave = np.sin(2 * np.pi * WAVE_FREQUENCY * time) * WAVE_STRETCH_Y / 2
    wave *= signal.windows.gaussian(len(time), WAVE_WIDTH)
    plt.plot(time, wave, linewidth=WAVE_LINEWIDTH, color=COLOR_3)

    # Format axes
    axes.set_axis_off()
    axes.set_aspect("equal")

    # Save figure
    figure.savefig(FILEPATH_LOGO_1)


def logo_square_small():

    # Figure and axes instance
    figure, axes = plt.subplots(1, figsize=SIZE_FIGURE_SMALL_INCHES)

    # Plot triangles
    for coordinates in zip(*item_coordinates_small()):
        x, y = coordinates
        if 2 - x == y:
            color = COLOR_1
        elif (x == 1) & (y == 1):
            color = COLOR_3
        else:
            color = COLOR_2
        axes.plot(
            *coordinates,
            "^",
            ms=SIZE_TRIANGLES * 0.8,
            clip_on=False,
            color=color
        )

    # Plot square braquets
    hspace = BRACKET_HORIZONTAL_MARGIN * 1.7
    wspace = BRACKET_VERTICAL_MARGIN * 1.7

    # Bottom left
    plt.plot(
        [-hspace, 0],
        [-wspace, -wspace],
        color=COLOR_3,
        linewidth=BRACKETS_LINEWIDTH,
    )

    # Bottom right
    plt.plot(
        [2 + hspace, 2],
        [-wspace, -wspace],
        color=COLOR_3,
        linewidth=BRACKETS_LINEWIDTH,
    )

    # Top left
    plt.plot(
        [-hspace, 0],
        [2 + wspace, 2 + wspace],
        color=COLOR_3,
        linewidth=BRACKETS_LINEWIDTH,
    )

    # Top right
    plt.plot(
        [2, 2 + hspace],
        [2 + wspace, 2 + wspace],
        color=COLOR_3,
        linewidth=BRACKETS_LINEWIDTH,
    )

    # Left
    plt.plot(
        [-hspace, -hspace],
        [-wspace, 2 + wspace],
        color=COLOR_3,
        linewidth=BRACKETS_LINEWIDTH,
    )

    # Right
    plt.plot(
        [2 + hspace, 2 + hspace],
        [-wspace, 2 + wspace],
        color=COLOR_3,
        linewidth=BRACKETS_LINEWIDTH,
    )

    # Wave
    time = np.linspace(-0.5, 0.5) * WAVE_STRETCH_X * 1.1
    wave = np.sin(2 * np.pi * WAVE_FREQUENCY / 2 * time) * WAVE_STRETCH_Y * 0.6
    wave *= signal.windows.gaussian(len(time), WAVE_WIDTH * 1.5)
    plt.plot(time, wave, linewidth=WAVE_LINEWIDTH, color=COLOR_3)

    # Text
    plt.text(
        4,
        1,
        "covseisnet",
        name="Ubuntu",
        va="center",
        ha="left",
        size=SIZE_TEXT * 1.5,
        color=COLOR_3,
        weight=500,
    )

    # Format axes
    axes.set_axis_off()
    axes.set_aspect("equal")

    # Save figur
    figure.savefig(FILEPATH_LOGO_2)


if __name__ == "__main__":
    logo_square()
    logo_square_small()
