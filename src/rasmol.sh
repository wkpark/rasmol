#!/bin/sh
RASMOLPATH=/usr/local/lib/rasmol
export RASMOLPATH

RASMOLPDBPATH=/data/brookhaven 
export RASMOLPDBPATH

rasmol.exe $*
