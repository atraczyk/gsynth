/*
*  Author: Andreas Traczyk <andreas.traczyk@savoirfairelinux.com>
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 3 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this program. If not, see <https://www.gnu.org/licenses/>.
*/
#pragma once

#include <iostream>
#include <iomanip>
#include <array>
#include <chrono>
#include <mutex>

#include <stdarg.h>

#undef min

#ifdef _DEBUG
#ifndef DBGOUT
#ifdef _WIN32
#define SEPERATOR '\\'
#else
#define SEPERATOR '/'
#endif // _WIN32
#ifndef _WIN32
#define LOGGER(m, ...)          consoleLog(LOG_FORMAT(m, ## __VA_ARGS__))
#define DBGOUT(m, ...)          consoleLog(m, ## __VA_ARGS__)
#else
#define FILENAME(X)              (strrchr(X, SEPERATOR) ? strrchr(X, SEPERATOR) + 1 : X)
#define STR(EXP)                #EXP
#define XSTR(X)                 STR(X)
//#define LINE                  
#define LOG_FORMAT(M, ...)      M, __VA_ARGS__
#define LOGGER(m, ...)          consoleLog(__FILE__, " :" XSTR(__LINE__), LOG_FORMAT(m, __VA_ARGS__))
#define DBGOUT(m, ...)          LOGGER(m, __VA_ARGS__)
#endif
#endif
#else
#define DBGOUT(m, ...)
#endif

static std::mutex logMutex;

inline void
printLog(std::ostream& s, const char* header, std::array<char, 8192>& buf, int len)
{
    std::lock_guard<std::mutex> lck(logMutex);

    // write timestamp
    using namespace std::chrono;
    using log_precision = milliseconds;
    constexpr auto den = log_precision::period::den;
    auto num = duration_cast<log_precision>(steady_clock::now().time_since_epoch()).count();
    s   << std::left << std::setw(20) << header 
        << "[" << std::setw(6) << num / den << "." << std::setw(3) << num % den << "]    "
        << std::right << std::fixed << buf.data();
    if ((size_t)len >= buf.size())
        s << "[[TRUNCATED]]";
    s << std::endl;
}

inline void
consoleLog(const char* filename, const char* lineno, char const *m, ...)
{
    std::array<char, 8192> buffer;
    va_list args;
    va_start(args, m);
    int ret = vsnprintf(buffer.data(), buffer.size(), m, args);
    va_end(args);

    if (ret < 0)
        return;

    printLog(   std::cout,
                std::string(FILENAME(filename)).append(lineno).c_str(),
                buffer,
                ret);
}