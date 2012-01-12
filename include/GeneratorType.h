#ifndef _ERHIC_MONTE_CARLO_GENERATOR_TYPE_H_
#define _ERHIC_MONTE_CARLO_GENERATOR_TYPE_H_

//
// GeneratorType.h
//
// Created by TB on 8/1/11.
// Copyright 2011 BNL. All rights reserved.
//

#include <iostream>
#include <string>

/**
 Describes types of eRHIC Monte Carlo generators and provides functions to
 determine which generator produced a particular output file.
*/
namespace erhic {
   
   namespace monte_carlo {
   
      class GeneratorType {
         
      public:
         
         enum {
            DJANGOH,
            MILOU,
            PEPSI,
            PYTHIA,
            RAPGAP,
            INVALID_GENERATOR
         };
         
         /**
          Attempt to determine which generator produced this output by reading
          the first line. Returns INVALID_GENERATOR if the generator is
          not valid, cannot be determined or if the stream is not good.
          After being called, the stream remains open after the first line.
          */
         static Int_t Determine(std::istream&);
         
         /**
          Attempt to determine which generator produced this output by reading
          the first line from the named file.
          Returns INVALID_GENERATOR if the generator is
          not valid, cannot be determined or if the file cannot be opened.
          */
         static Int_t Determine(const std::string&);
      };
      
   } // namespace monte_carlo
   
} // namespace erhic

#endif
