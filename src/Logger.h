#pragma once

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

class Logger {
public:
    Logger();  // Constructor to initialize logger
    ~Logger(); // Destructor

    // Method to initialize the logger with a specified log level
    void initLogger(const std::string& log_level);

    // Methods to log messages at different levels
    void trace(const std::string &message);
    void debug(const std::string& message);
    void info(const std::string& message);
    void warn(const std::string& message);
    void error(const std::string& message);
    void critical(const std::string& message);




private:
    std::shared_ptr<spdlog::logger> logger;  // Logger instance
};

