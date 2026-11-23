//
// Created by xos81802 on 26/07/2026.
//

#include "Mask.h"


// Mask - read in the coordinates of the particle (+ hydration)
// mapped waters mark the location of potential surface
// build grid - particle is centered.
// grid should extend additional 7.2 Angstroms away from protein
// remove grid points that contain protein
// empty grid points represent putative waters

// calculate cross-term with particle to determine weight
