#!/usr/bin/env python3
"""Make a two-slice CRAM container from the checked-in single-slice fixture."""
import gzip
from pathlib import Path
import struct
import zlib


def read_itf8(data, position):
    first = data[position]
    size = next((n for n in range(1, 5) if first < (0x100 - (0x80 >> (n - 1)))), 5)
    if size == 5:
        value = ((first & 15) << 28) | (data[position + 1] << 20) | (data[position + 2] << 12) | (data[position + 3] << 4) | (data[position + 4] & 15)
    else:
        value = first & (0xff >> size)
        for byte in data[position + 1:position + size]:
            value = (value << 8) | byte
    return value, position + size


def itf8(value):
    value &= 0xffffffff
    for size in range(1, 5):
        if value < (1 << (7 * size)):
            encoded = bytearray(value.to_bytes(size, 'big'))
            encoded[0] |= (0xff << (9 - size)) & 0xff
            return bytes(encoded)
    return bytes([0xf0 | (value >> 28), (value >> 20) & 255, (value >> 12) & 255, (value >> 4) & 255, value & 15])


def read_ltf8(data, position):
    first = data[position]
    leading = 0
    while leading < 8 and first & (0x80 >> leading):
        leading += 1
    value = first & (0x7f >> leading)
    for byte in data[position + 1:position + leading + 1]:
        value = (value << 8) | byte
    return value, position + leading + 1


def ltf8(value):
    if value < 128:
        return bytes([value])
    raise ValueError('This fixture only needs small record counters')


def container(data, offset):
    length = struct.unpack_from('<i', data, offset)[0]
    position = offset + 4
    fields = []
    for _ in range(4):
        value, position = read_itf8(data, position)
        fields.append(value)
    for _ in range(2):
        value, position = read_ltf8(data, position)
        fields.append(value)
    blocks, position = read_itf8(data, position)
    count, position = read_itf8(data, position)
    landmarks = []
    for _ in range(count):
        value, position = read_itf8(data, position)
        landmarks.append(value)
    assert struct.unpack_from('<I', data, position)[0] == zlib.crc32(data[offset:position])
    position += 4
    return fields, blocks, landmarks, data[position:position + length]


def recounter_slice(data, counter):
    method, content_type = data[:2]
    assert content_type == 2
    content_id, position = read_itf8(data, 2)
    compressed_size, position = read_itf8(data, position)
    raw_size, position = read_itf8(data, position)
    end = position + compressed_size
    assert struct.unpack_from('<I', data, end)[0] == zlib.crc32(data[:end])
    payload = data[position:end]
    if method == 1:
        payload = gzip.decompress(payload)
    else:
        assert method == 0
    assert len(payload) == raw_size
    position = 0
    for _ in range(4):
        _, position = read_itf8(payload, position)
    _, counter_end = read_ltf8(payload, position)
    payload = payload[:position] + ltf8(counter) + payload[counter_end:]
    block = bytes([0, content_type]) + itf8(content_id) + itf8(len(payload)) * 2 + payload
    return block + struct.pack('<I', zlib.crc32(block)) + data[end + 4:]


def main():
    directory = Path(__file__).resolve().parent
    source = directory / 'single_ref_slices.cram'
    data = source.read_bytes()
    assert data[:5] == b'CRAM\x03' and data[5] in (0, 1)
    entries = [list(map(int, line.split('\t'))) for line in gzip.decompress(source.with_suffix('.cram.crai').read_bytes()).decode().splitlines()]
    offset = next(entry[3] for entry in entries if entry[0] == 0)
    unplaced = next(entry[3] for entry in entries if entry[0] == -1)
    fields, blocks, landmarks, payload = container(data, offset)
    assert len(landmarks) == 1 and fields[0] == 0 and fields[3] == 2 and fields[4] == 0
    header_block = payload[:landmarks[0]]
    first_slice = payload[landmarks[0]:]
    second_slice = recounter_slice(first_slice, 2)
    payload = header_block + first_slice + second_slice
    fields[3] *= 2
    fields[5] *= 2
    landmarks.append(landmarks[0] + len(first_slice))
    header = struct.pack('<i', len(payload)) + b''.join(itf8(n) for n in fields[:4])
    header += b''.join(ltf8(n) for n in fields[4:])
    header += itf8(1 + 2 * (blocks - 1)) + itf8(2) + b''.join(itf8(n) for n in landmarks)
    output = data[:offset] + header + struct.pack('<I', zlib.crc32(header)) + payload + data[unplaced:]
    (directory / 'same_reference_multi_slice.cram').write_bytes(output)


if __name__ == '__main__':
    main()
