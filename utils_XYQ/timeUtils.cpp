

#include <chrono>
#include <iostream>
#include <string>
#include <vector>

#include "timeUtils.h"

using namespace std;


namespace timeUtils{

    string get_ymd_date_str( ymd_date tarDate ){
        string ymd_date_str = to_string( tarDate.year ) + "/" +
            to_string( tarDate.month ) + "/" +
             to_string( tarDate.day );
        return ymd_date_str;
    }

    void print_ymd_date( ymd_date tarDate ){
        cout << timeUtils::get_ymd_date_str( tarDate ) << endl;
    }

}


chrono::system_clock::duration getDuration(){

    // Get current time
    auto now = chrono::system_clock::now();
    // Convert the current time to time since epoch
    auto duration = now.time_since_epoch();

    return duration;

}


double getSecSinceUnix(){
    
    // Convert the current time to time since epoch
    auto duration = getDuration();
    // Convert duration to seconds
    auto secs
        = chrono::duration_cast<chrono::seconds>(duration).count();
    return (double) secs;

}


double getMSecSinceUnix(){
    
    // Convert the current time to time since epoch
    auto duration = getDuration();
    // Convert duration to milli-seconds
    auto msecs
        = chrono::duration_cast<chrono::milliseconds>(duration).count();
    return (double) msecs;

}