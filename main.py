#! /usr/bin/env python

import os

from redshift import cluster

def parse_config():
    """Parse the config file, making a dictionary of parameters

    :return: dictionary, with keys of parameter names, and values of the
    desired values for those parameters.
    """
    global params  # params needs to be available elsewhere
    params = dict()
    config_file = open("config.txt")
    for line in config_file:
        if not line.startswith("#") and not line.isspace():
            params[line.split()[0]] = line.split()[2]
    config_file.close()
    return params


def main():
    # parse the config file
    params = parse_config()

    catalogs = os.listdir(params["catalog_directory"])
    for cat in catalogs:
        filepath = params["catalog_directory"] + cat
        cl = cluster.Cluster(filepath)




if __name__ == "__main__":
    main()


# TODO: What is the format of the catalogs I will get?