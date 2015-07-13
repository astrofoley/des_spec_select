#!/bin/sh
sed 's/,//g;s/"//g;s/\[//g;s/]//g;s/}//g' out_temp > out_temp3
