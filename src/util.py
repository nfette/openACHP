def optimize_result_repr_to_dict(s):
    """Warning: dangerous function. Uses eval()."""
    import re
    from numpy import array, nan, inf
    ss = re.split(r'\s+(\w+): ', s)

    i = 1
    d = dict()
    while i < len(ss):
        name = ss[i]
        i += 1
        value = ss[i]
        i += 1
        d[name] = eval(value)
    return d
