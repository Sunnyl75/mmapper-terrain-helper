"""Read MMapper room geometry from XML or supported MM2 room records."""
from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
import re
import struct
import unicodedata
import xml.etree.ElementTree as ET
import zipfile
import zlib

TERRAIN=("UNDEFINED","INDOORS","CITY","FIELD","FOREST","HILLS","MOUNTAINS",
         "SHALLOW","WATER","RAPIDS","UNDERWATER","ROAD","BRUSH","TUNNEL","CAVERN")
OUTDOOR=frozenset(("FIELD","FOREST","HILLS","MOUNTAINS","SHALLOW","WATER","RAPIDS","ROAD","BRUSH"))
CARDINAL=("north","south","east","west")
DIRS=CARDINAL+("up","down","unknown")
STEPS=((0,1),(0,-1),(1,0),(-1,0))
END=0xffffffff


@dataclass(frozen=True,slots=True)
class Room:
    id: int
    server_id: int
    xyz: tuple
    terrain: str
    name: str
    description: str
    flags: int
    exits: dict

    @property
    def outdoor(self):
        return self.xyz[2]==0 and self.terrain in OUTDOOR


def normalize(text):
    return re.sub(r"\s+"," ",unicodedata.normalize("NFKC",text).casefold()).strip()


def fingerprint(room):
    desc=normalize(room.description)
    return (normalize(room.name),sha256(desc.encode()).digest()) if len(desc)>=56 and room.name else None


def xml_rooms(stream):
    events=iter(ET.iterparse(stream,events=("start","end")))
    _,root=next(events)
    if root.tag!="map" or root.get("type")!="mmapper2xml":
        raise ValueError("Expected an exported MMapper XML map")
    for event,e in events:
        if event!="end" or e.tag!="room":
            continue
        c=e.find("coord")
        if c is None:
            raise ValueError(f"Room {e.get('id')} has no coordinate")
        flags=0
        exits={}
        for x in e.findall("exit"):
            key=x.get("dir")
            exits[key]=tuple(int(t.text) for t in x.findall("to") if t.text)
            if key in CARDINAL and any(f.text=="ROAD" for f in x.findall("exitflag")):
                flags|=1<<CARDINAL.index(key)
        yield Room(int(e.get("id")),int(e.get("server_id") or 0),
                   tuple(int(c.get(k)) for k in "xyz"),e.findtext("terrain",default="UNDEFINED"),
                   e.get("name",""),e.findtext("description",default=""),flags,exits)
        root.remove(e)


class QtRoomStream:
    def __init__(self,data):
        self.data=data
        self.pos=0

    def unpack(self,code):
        fmt=struct.Struct(">"+code)
        value=fmt.unpack_from(self.data,self.pos)
        self.pos+=fmt.size
        return value[0] if len(value)==1 else value

    def string(self,keep=True):
        n=self.unpack("I")
        if n==END:
            return ""
        if n%2 or n>len(self.data)-self.pos:
            raise ValueError("Invalid QString in MM2 map")
        p=self.pos
        self.pos+=n
        return self.data[p:p+n].decode("utf-16-be",errors="replace") if keep else ""

    def links(self):
        links=[]
        while True:
            i=self.unpack("I")
            if i==END:
                return tuple(links)
            links.append(i)


def binary_rooms(path):
    raw=path.read_bytes()
    try:
        magic,version,length=struct.unpack_from(">III",raw)
        if magic!=0xffb2af01 or not 36<=version<=42:
            raise ValueError("Supported MM2 schemas are 36–42. Export other versions to MMapper XML.")
        data=zlib.decompress(raw[12:])
        if len(data)!=length:
            raise ValueError("MM2 decompressed length does not match its header")
        f=QtRoomStream(data)
        count,_,_,_,_=f.unpack("IIiii")
        for _ in range(count):
            if version>=42:
                f.string(False)
            name,desc=f.string(),f.string()
            f.string(False) # contents
            id_=f.unpack("I")
            server_id=f.unpack("I") if version>=40 else 0
            f.string(False) # note
            terrain=f.unpack("B")
            f.unpack("5B2I") # light, align, portable, ridable, sundeath, mob/load flags
            if version<39:
                f.unpack("B")
            xyz=f.unpack("iii")
            exits={}
            flags=0
            for index,direction in enumerate(DIRS):
                bits,door=f.unpack("HH")
                f.string(False)
                if version<38:
                    f.links()
                exits[direction]=f.links()
                if index<4 and bits&4:
                    flags|=1<<index
            kind=TERRAIN[terrain] if terrain<len(TERRAIN) else "INDOORS" if terrain==15 and version<41 else "UNDEFINED"
            yield Room(id_,server_id,xyz,kind,name,desc,flags,exits)
    except (IndexError,struct.error,zlib.error) as error:
        raise ValueError(f"Invalid or unsupported room data in {path.name}") from error


def read_map(path):
    path=Path(path)
    if zipfile.is_zipfile(path):
        with zipfile.ZipFile(path) as archive:
            maps=[s for s in archive.namelist() if s.endswith((".mm2xml",".xml"))]
            if len(maps)!=1:
                raise ValueError("A grid archive must contain one MMapper XML file")
            with archive.open(maps[0]) as input_:
                return list(xml_rooms(input_))
    with path.open("rb") as input_:
        if input_.read(4)==b"\xff\xb2\xaf\x01":
            return list(binary_rooms(path))
        input_.seek(0)
        return list(xml_rooms(input_))
