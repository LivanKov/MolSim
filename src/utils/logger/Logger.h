#pragma once

#include <memory>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>
#include <string>

/**
 * @class Logger.
 * @brief A (singleton) Logger class for logging operations though out the
 * application. This class provides various logging methods for any specific
 * log-level (trace, debug, info, warn, error, off). Using the singleton pattern
 * to ensure that only one instance of the logger is created. The class uses the
 * spdlog library.
 */
class Logger {
public:
  /**
   * @brief To retrieve the logger instance.
   * @param Sets the default log-level to 'info'.
   * @return A reference to the logger instance.
   */
  static Logger &getInstance(const std::string &log_level = "info");

  /**
   * @brief Logs message with respective log-level.
   * @param Logging message.
   */
  void trace(const std::string &message);
  void debug(const std::string &message);
  void info(const std::string &message);
  void warn(const std::string &message);
  void error(const std::string &message);

private:
  /**
   * @brief Creates Logger instance with specified log-level.
   * @param Sets the initial log-level.
   *
   * This constructor is private to enforce singleton pattern.
   * Prevents external instantiation.
   */
  Logger(const std::string &log_level);

  /**
   * @brief Destructor for the Logger class.
   */
  ~Logger();

  std::shared_ptr<spdlog::logger> logger;

  // No instance of Logger class can be copied
  Logger(const Logger &) = delete;

  // Prevents assigning instance of Logger class to another instance
  Logger &operator=(const Logger &) = delete;
};
