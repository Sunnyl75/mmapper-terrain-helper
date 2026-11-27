"""Write and check terrain-only MMapper schema 42 backgrounds."""
import struct
import zlib

from room_io import DIRS, END, OUTDOOR, STEPS, TERRAIN, binary_rooms

MAGIC = 0xFFB2AF01
SCHEMA = 42
ROOM = struct.Struct(">7I6B2I3i")
EXIT = struct.Struct(">HHI")
UINT = struct.Struct(">I")


def write_mm2(path, grid, flags, bounds):
    west, south, east, north = map(int, bounds)
    width, height = east-west, north-south
    if grid.shape != (height, width) or flags.shape != grid.shape:
        raise ValueError("Background grid has unexpected dimensions")
    buffer = bytearray(struct.pack(">IIiii", width*height, 0, west, north-1, 0))
    for y in range(north-1, south-1, -1):
        for x in range(west, east):
            key = (north-1-y)*width+x-west
            cell = y-south, x-west
            # Four empty QStrings, synthetic ID, no server ID, empty note.
            buffer.extend(ROOM.pack(0,0,0,0,key,0,0,int(grid[cell]),0,0,0,0,0,0,0,x,y,0))
            mask = int(flags[cell])
            for i, (dx,dy) in enumerate(STEPS):
                exists = west<=x+dx<east and south<=y+dy<north
                buffer.extend(EXIT.pack(int(exists) | (4 if mask & (1<<i) else 0),0,0))
                if exists:
                    buffer.extend(UINT.pack(key+dx-dy*width))
                buffer.extend(UINT.pack(END))
            for _ in DIRS[len(STEPS):]:
                buffer.extend(EXIT.pack(0,0,0))
                buffer.extend(UINT.pack(END))
    path.write_bytes(struct.pack(">III", MAGIC, SCHEMA, len(buffer))+zlib.compress(buffer,9))


def verify_mm2(path, grid, flags, bounds):
    west, south, east, north = map(int, bounds)
    width, height = east-west, north-south
    size = 0
    for i, room in enumerate(binary_rooms(path)):
        x,y,z = room.xyz
        if room.id!=i or z or x!=west+i%width or y!=north-1-i//width:
            raise ValueError("MM2 coordinate or ordering mismatch")
        cell = y-south,x-west
        if room.terrain not in OUTDOOR or room.terrain!=TERRAIN[grid[cell]] or room.flags!=flags[cell]:
            raise ValueError(f"MM2 terrain or road mask mismatch at {(x,y)}")
        if room.name or room.description or room.server_id:
            raise ValueError("MM2 contains foreground gameplay metadata")
        for j,direction in enumerate(DIRS):
            expected=()
            if j<len(STEPS):
                dx,dy=STEPS[j]
                if west<=x+dx<east and south<=y+dy<north:
                    expected=(i+dx-dy*width,)
            if room.exits[direction]!=expected:
                raise ValueError("MM2 contains a non-grid exit")
        size=i+1
    if size!=width*height:
        raise ValueError("MM2 has incomplete background coverage")
    return {"schema":SCHEMA,"bytes":path.stat().st_size,"verified_cells":size}
