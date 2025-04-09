#pragma once

#include <chrono>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

namespace timeUtils{

    // This structure is a year/month/day format date.
    struct ymd_date { 
        int year = 0;
        int month = 0;
        int day = 0;
    };

    /* Obtain string of the date. */
    string get_ymd_date_str( ymd_date tarDate );

    /* Print the date. */
    void print_ymd_date( ymd_date tarDate );

}



chrono::system_clock::duration getDuration();

/* 
    Obtain the number of seconds since Epoch.
*/
double getSecSinceUnix();

/* 
    Obtain the number of milliseconds since Epoch.
*/
double getMSecSinceUnix();


