#!/usr/bin/env bash

make dirs
tar -xzf interpolator-data interpolator-data.tar.gz
./eic $@
