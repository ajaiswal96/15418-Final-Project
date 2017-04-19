#include <stdio.h>
#include <stdlib.h>

#include "serial.h"

int main (int argc, char** argv) {
    printf("starting main\n");
    serial(1);
    printf("done with main\n");
}
