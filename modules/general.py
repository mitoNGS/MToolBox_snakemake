#!/usr/bin/env python

def memory_usage_resource():
    """
    Function to get memory usage (in MB)
    Source: http://fa.bianp.net/blog/2013/different-ways-to-get-memory-consumption-or-lessons-learned-from-memory_profiler/
    """
    import resource
    rusage_denom = 1024.
    if sys.platform == 'darwin':
        # ... it seems that in OSX the output is different units ...
        rusage_denom = rusage_denom * rusage_denom
    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / rusage_denom
    return mem

def s_encoding(s):
    if type(s) == bytes:
        return s.decode("utf-8")
    elif type(s) == str:
        return s

def get_SAM_header(samfile):
    # is the file compressed?
    if samfile.endswith("gz"):
        samhandle = gzip.GzipFile(samfile, mode = 'r')
    elif samfile.endswith("bz2"):
        samhandle = bz2.BZ2File(samfile, mode = 'r')
    else:
        samhandle = open(samfile, 'r')
    comment_count = 0
    header_lines = []
    l = s_encoding(samhandle.readline())
    print(l)
    while l[0] == "@":
        header_lines.append(l)
        comment_count += 1
        l = s_encoding(samhandle.readline())
    return header_lines, comment_count
