##
# Cloth
#
# @file
# @version 0.1

all:
	cc -o cloth main.c -lraylib -lGL -lm -lpthread -ldl -lrt -lX11 && ./cloth
# end
