#include "stplugin.h"
#include "honesttest.h"

STDLL stata_call(int argc, char *argv[])
{
    SF_display("Hello World\n");
    return(0);
}
