#ifndef MAIN_H
#define MAIN_H

#include <errno.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <sys/wait.h>
#include <sys/types.h>

#include <iostream>
                              
#include <Manzoni.h>
#include <utils/signals.h>
#include <utils/Logger.h>
#include <utils/Arguments.h>  
#include <utils/TypeConversions.h>

// Messages
void startMessage();
void helpMessage(bool isModuleSelected = true); 
void compilationOptions();      
void version();    

// Helpers
void signalsDefinition();

#endif
