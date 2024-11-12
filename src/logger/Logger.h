#pragma once

#include <memory>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>
#include <string>

class Logger {
public:
  static Logger &getInstance(const std::string &log_level = "info");

  // Logging methods
  void trace(const std::string &message);
  void debug(const std::string &message);
  void info(const std::string &message);
  void warn(const std::string &message);
  void error(const std::string &message);
  void critical(const std::string &message);

private:
  // Prevents external instantiation
  Logger(const std::string &log_level);


  ~Logger();

  std::shared_ptr<spdlog::logger> logger;

  // No instance of Logger class can be copied
  Logger(const Logger &) = delete;
  // Prevents assigning instance of Logger class to another instance
  Logger &operator=(const Logger &) = delete;
};
