#include "Logger.h"

Logger::Logger(const std::string &log_level) {
  // Create a color console sink for coloured logs
  auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
  console_sink->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%l] %v");
  console_sink->set_pattern("%+");

  // Initialize the logger with the sink
  logger = std::make_shared<spdlog::logger>("console", console_sink);

  // Set the log level based on user input
  if (log_level == "trace") {
    logger->set_level(spdlog::level::trace);
  } else if (log_level == "debug") {
    logger->set_level(spdlog::level::debug);
  } else if (log_level == "info") {
    logger->set_level(spdlog::level::info);
  } else if (log_level == "warn") {
    logger->set_level(spdlog::level::warn);
  } else if (log_level == "error") {
    logger->set_level(spdlog::level::err);
  } else if (log_level == "off") {
    logger->set_level(spdlog::level::off);
  } else {
    logger->set_level(spdlog::level::info);
  }
}

Logger::~Logger() {}

Logger &Logger::getInstance(const std::string &log_level) {
  static Logger instance(log_level); // Ensures to be initialized only once
  return instance;
}

void Logger::trace(const std::string &message) { logger->trace(message); }

void Logger::debug(const std::string &message) { logger->debug(message); }

void Logger::info(const std::string &message) { logger->info(message); }

void Logger::warn(const std::string &message) { logger->warn(message); }

void Logger::error(const std::string &message) { logger->error(message); }
