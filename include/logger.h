#ifndef NEGATIVE_TRIANGLE_COUNTING_LOGGER_H
#define NEGATIVE_TRIANGLE_COUNTING_LOGGER_H

#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/console.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>


inline void init_logging(const std::string* file_name) {
    namespace logging = boost::log;
    namespace keywords = boost::log::keywords;

    logging::add_common_attributes();

    if(file_name == nullptr) {
        logging::add_console_log(std::cout);
    } else {
        logging::add_file_log(
            keywords::file_name = *file_name,
            keywords::auto_flush = true
        );
    }
}

#endif //NEGATIVE_TRIANGLE_COUNTING_LOGGER_H