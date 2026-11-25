"""Generate the small procedural OBJ used for command-line packaging tests."""
from pathlib import Path
import sys


def grid():
    side = 12
    lines = [f'v {x / side} {y / side} 0'
             for y in range(side + 1) for x in range(side + 1)]
    for y in range(side):
        for x in range(side):
            first = y * (side + 1) + x + 1
            below = first + side + 1
            lines.extend((f'f {first} {first + 1} {below + 1}',
                          f'f {first} {below + 1} {below}'))
    return '\n'.join(lines) + '\n'


if __name__ == '__main__':
    Path(sys.argv[1]).write_text(grid(), encoding='ascii')
