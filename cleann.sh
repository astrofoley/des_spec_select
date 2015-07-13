#!/bin/sh
sed 's/,//g;s/"//g;s/\[//g;s/]//g;s/}//g' out_temp2 > out_temp3
