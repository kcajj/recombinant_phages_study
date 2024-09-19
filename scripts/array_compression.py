def compress_array(array):
    compressed_array = []
    previous = array[0]
    c = 0
    for i in range(len(array)):
        if array[i] == previous:
            c += 1
        else:
            compressed_array.append((c, int(previous)))
            c = 1
        previous = array[i]
    compressed_array.append((c, int(previous)))
    return compressed_array

def decompress_array(compressed_array):
    array = []
    for c, v in compressed_array:
        array += [int(v)] * c
    return array