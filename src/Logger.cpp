/* #include "Logger.h"

Logger::Logger() {
    // Default constructor, the logger is not initialized yet.
}

Logger::~Logger() {
    // Destructor, cleanup if needed.
}

void Logger::initLogger(const std::string& log_level) {
    // Create a color console sink for pretty logs
    auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
    // console_sink->set_pattern("%^[%T] %n: %v%$"); // Log format with timestamp, logger name, and message
    console_sink->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%l] %v");

    // console_sink->set_pattern("%+");

    // Initialize the logger with the sink
    logger = std::make_shared<spdlog::logger>("console", console_sink);

    // Set the log level based on user input
    if (log_level == "trace") {
        logger->set_level(spdlog::level::trace); // Set to trace level
    } else if (log_level == "debug") {
        logger->set_level(spdlog::level::debug); // Set to debug level
    } else if (log_level == "info") {
        logger->set_level(spdlog::level::info); // Set to info level
    } else if (log_level == "warn") {
        logger->set_level(spdlog::level::warn); // Set to warn level
    } else if (log_level == "error") {
        logger->set_level(spdlog::level::err); // Set to error level
    } else if (log_level == "critical") {
        logger->set_level(spdlog::level::critical); // Set to critical level
    } else if (log_level == "off") {
        logger->set_level(spdlog::level::off); // Set to off level
    }
    else {
        logger->set_level(spdlog::level::info); // Default to info level
    }
}


void Logger::trace(const std::string& message) {
    logger->trace(message);
}

void Logger::debug(const std::string& message) {
    logger->debug(message);
}

void Logger::info(const std::string& message) {
    logger->info(message);
}

void Logger::warn(const std::string& message) {
    logger->warn(message);
}

void Logger::error(const std::string& message) {
    logger->error(message);
}

void Logger::critical(const std::string& message) {
    logger->critical(message);
} */

#include "Logger.h"

Logger::Logger(const std::string& log_level) {
    // Create a color console sink for pretty logs
    auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
    console_sink->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%l] %v");
    console_sink->set_pattern("%+");

    // Initialize the logger with the sink
    logger = std::make_shared<spdlog::logger>("console", console_sink);

    // Set the log level based on user input
    if (log_level == "trace") {
        logger->set_level(spdlog::level::trace); // Set to trace level
    } else if (log_level == "debug") {
        logger->set_level(spdlog::level::debug); // Set to debug level
    } else if (log_level == "info") {
        logger->set_level(spdlog::level::info); // Set to info level
    } else if (log_level == "warn") {
        logger->set_level(spdlog::level::warn); // Set to warn level
    } else if (log_level == "error") {
        logger->set_level(spdlog::level::err); // Set to error level
    } else if (log_level == "critical") {
        logger->set_level(spdlog::level::critical); // Set to critical level
    } else if (log_level == "off") {
        logger->set_level(spdlog::level::off); // Set to off level
    } else {
        logger->set_level(spdlog::level::info); // Default to info level
    }
}

Logger::~Logger() {
    // Optional: Clean-up if necessary, but spdlog will manage memory for us
}

Logger& Logger::getInstance(const std::string& log_level) {
    static Logger instance(log_level);  // Guaranteed to be initialized only once
    return instance;
}

void Logger::trace(const std::string& message) {
    logger->trace(message);
}

void Logger::debug(const std::string& message) {
    logger->debug(message);
}

void Logger::info(const std::string& message) {
    logger->info(message);
}

void Logger::warn(const std::string& message) {
    logger->warn(message);
}

void Logger::error(const std::string& message) {
    logger->error(message);
}

void Logger::critical(const std::string& message) {
    logger->critical(message);
}

