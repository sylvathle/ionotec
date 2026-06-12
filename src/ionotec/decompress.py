import sys
import os

# ---------------------------------------------------------------------------
# LZW decompression (Unix compress format)
# ---------------------------------------------------------------------------
 
MAGIC_0    = 0x1F
MAGIC_1    = 0x9D
BIT_MASK   = 0x1F   # lower 5 bits of header[2] → max code-width
BLOCK_FLAG = 0x80   # bit 7 of header[2]  → block-compress mode
 
 
def decompress_z(data: bytes) -> bytes:
    """Decompress a bytes object containing a .Z (Unix LZW) payload.
 
    Raises ValueError on bad input.
    Returns the decompressed data as bytes.
    """
    # ------------------------------------------------------------------ header
    if len(data) < 3 or data[0] != MAGIC_0 or data[1] != MAGIC_1:
        raise ValueError("Not a valid .Z file (bad magic bytes)")
 
    flags      = data[2]
    if flags & 0x60:
        raise ValueError("Invalid header flags (reserved bits set)")
 
    max_bits   = flags & BIT_MASK
    if not (9 <= max_bits <= 16):
        raise ValueError(f"Unsupported max_bits: {max_bits} (must be 9–16)")
    if max_bits == 9:
        max_bits = 10          # 9 is treated as 10 by the original tool
 
    block_mode = bool(flags & BLOCK_FLAG)
 
    if len(data) == 3:
        return b""             # empty payload is valid
 
    # ------------------------------------------------------ initial state
    bits  = 9
    mask  = 0x1FF              # (1 << bits) - 1
    end   = 256 if block_mode else 255   # highest code in table so far
 
    # String table stored as parallel prefix/suffix arrays (fast & classic)
    prefix = [0] * 65536
    suffix = [0] * 65536
 
    # -------------------------------------------------------- bootstrap
    # Read the very first code (must be a literal 0–255) without building
    # a table entry yet.
    buf   = data[3] | (data[4] << 8)
    final = prev = buf & mask
    buf >>= bits
    left  = 16 - bits          # bits still valid in buf
 
    if prev > 255:
        raise ValueError("First code is not a literal")
 
    put  = [final]             # output accumulator
    mark = 3                   # byte offset where current code-width block started
    nxt  = 5                   # next byte to consume from data[]
 
    # ---------------------------------------------------- main decode loop
    while nxt < len(data):
 
        # --- widen code if the table is about to outgrow the current mask ---
        if end >= mask and bits < max_bits:
            # Flush to the next 8*bits-byte boundary from mark
            # (legacy behaviour inherited from a VAX implementation)
            rem = (nxt - mark) % bits
            if rem:
                rem = bits - rem
                if rem >= len(data) - nxt:
                    break
                nxt += rem
 
            buf  = 0
            left = 0
            mark = nxt         # new reference point for the next flush
 
            bits += 1
            mask  = (mask << 1) | 1
 
        # --- read one code of `bits` bits ---
        buf  += data[nxt] << left
        nxt  += 1
        left += 8
        if left < bits:
            if nxt == len(data):
                raise ValueError("Stream ended mid-code")
            buf  += data[nxt] << left
            nxt  += 1
            left += 8
 
        code  = buf & mask
        buf >>= bits
        left -= bits
 
        # --- handle CLEAR code (block mode only) ---
        if block_mode and code == 256:
            # Flush to next 8*bits-byte boundary
            rem = (nxt - mark) % bits
            if rem:
                rem = bits - rem
                if rem > len(data) - nxt:
                    break
                nxt += rem
            buf  = 0
            left = 0
            mark = nxt
 
            bits  = 9
            mask  = 0x1FF
            end   = 255        # empty table (literals only)
            continue
 
        # --- decode the code to a string ---
        temp  = code
        stack = []
 
        if code > end:
            # Code not yet in table: must be the "ω+ω[0]" special case
            if code != end + 1 or prev > end:
                raise ValueError(f"Invalid LZW code {code} (end={end})")
            stack.append(final)
            code = prev
 
        while code >= 256:
            stack.append(suffix[code])
            code = prefix[code]
 
        stack.append(code)
        final = code           # first character of decoded string
 
        # --- add new entry to the string table ---
        if end < mask:
            end += 1
            prefix[end] = prev
            suffix[end] = final
 
        prev = temp
 
        # --- emit decoded string (stack is reversed) ---
        put += stack[::-1]
 
    return bytes(bytearray(put))
 