

def read_catalog(file_path):
    """
    Read a catalog at the given location, and turn it into a cluster object.

    The catalog should be formatted in the way that all of them are.
    TODO: what is this format?

    :param file_path: complete file path of the catalog with photometric
                    data of the cluster.
    :return: cluster object with data from the catalog
    """
    cat_file = open(file_path)


def redshift(file_path):
    """

    :param file_path: complete path of the catalog containing the photometric
                    data of the cluster.
    :return: data object, with the redshift and errors
    """

    # read the data in the catalog into a cluster object.
    cluster = read_catalog(file_path)