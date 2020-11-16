#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import yaml
from collections import OrderedDict

inconfig = snakemake.config


def od2dict(d):
    if isinstance(d, dict):
        tmp = dict()
        for k, v in d.items():
            tmp[k] = od2dict(v)
        return tmp
    return d


config = od2dict(inconfig)

with open(snakemake.output[0], "w") as fh:
    fh.write(yaml.dump(config))
