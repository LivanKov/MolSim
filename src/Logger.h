/* #pragma once

#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

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
}; */

#pragma once

#include <memory>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>
#include <string>

class Logger {
public:
  // Public static method to get the logger instance
  static Logger &getInstance(const std::string &log_level = "info");

  // Logging methods
  void trace(const std::string &message);
  void debug(const std::string &message);
  void info(const std::string &message);
  void warn(const std::string &message);
  void error(const std::string &message);
  void critical(const std::string &message);

private:
  // Private constructor to prevent external instantiation
  Logger(const std::string &log_level);

  // Private destructor (optional, if you want more control)
  ~Logger();

  // Shared pointer to the actual spdlog logger
  std::shared_ptr<spdlog::logger> logger;

  // Delete copy constructor and assignment operator
  Logger(const Logger &) = delete;
  Logger &operator=(const Logger &) = delete;
};
