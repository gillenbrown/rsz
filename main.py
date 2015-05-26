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

def start_output():
    """Start the output file, and put the header in the file.

    :return: open file object pointing to the output file.
    """
    output_file = open(params["results_file"], "w")
    header = '# {:15s}{:6s}{:6s}{:6s}'.format("name", "z", "ezu", "ezd")
    output_file.write(header)
    return output_file

def add_to_catalog(output_file, cluster):
    """

    :param cluster:
    :return:
    """
    line = '\n{:17s}{:<6f}{:<6f}{:<6f}'.format(cluster.name, cluster.z.value,
                                          cluster.z.upper_error,
                                          cluster.z.lower_error)
    output_file.write(line)


def main():
    # parse the config file
    params = parse_config()

    catalogs = os.listdir(params["catalog_directory"])
    # weed out things that aren't catalogs
    catalogs = [cat for cat in catalogs if cat.endswith(params["extension"])]

    # start the output file
    output_file = start_output()

    for cat in catalogs:
        filepath = params["catalog_directory"] + cat
        cl = cluster.Cluster(filepath)

        # do the fitting procedure
        cl.fit_z(params)

        # make the output catalog
        add_to_catalog(output_file, cl)

    output_file.close()




if __name__ == "__main__":
    main()


# TODO: What is the format of the catalogs I will get?