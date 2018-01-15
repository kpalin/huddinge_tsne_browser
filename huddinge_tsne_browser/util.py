def memory_usage_psutil():
    # return the memory usage in MB
    import psutil
    import os
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / float(2**20)
    return mem


def memory_usage_resource():
    import resource
    import sys
    rusage_denom = 1024.
    if sys.platform == 'darwin':
        # ... it seems that in OSX the output is different units ...
        rusage_denom = rusage_denom * rusage_denom
    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / rusage_denom
    return mem


def memory_usage():
    "Report memory usage of current process"
    import logging as log
    try:
        mem = memory_usage_psutil()
    except ImportError as e:
        log.info(
            "No psutil found. Reporting memory usage from resource module: %s"
            % (str(e)))
        mem = memory_usage_resource()

    return mem
