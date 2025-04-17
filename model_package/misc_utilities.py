
def str2bool(v):
    '''Function for converting string to Boolean. Borrowed from: https://stackoverflow.com/questions/15008758/parsing-boolean-values-with-argparse

    :param str/bool v: A string or boolean indicating a True or False value

    :returns: True or False
    '''

    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')