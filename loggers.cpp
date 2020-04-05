#include "loggers.h"
#include <string>

namespace Loggers
{

Logger::Logger()
{
    level = LogLevel::Debug;
}

bool Logger::should_log (const LogLevel& message_level) const
{
    return (message_level >= level);
}

void StdoutLogger::log (const LogLevel& message_level, const std::string& message) const
{
    printf("%s", message.c_str());
}

CombinedLogger::CombinedLogger()
{
    loggers = std::make_unique<std::vector<Logger>>();
}

void CombinedLogger::log (const LogLevel& message_level, const std::string& message) const
{
    for (auto& it : *loggers)
    {
        it.log(message_level, message);
    }
}

std::unique_ptr<Logger> LoggerFactory::make_stdout_logger(const LogLevel& level)
{
    std::unique_ptr<Logger> logger = std::make_unique<StdoutLogger>();
    logger->level = level;
    return logger;
}

NullLogger* LoggerFactory::null_logger_instance = new NullLogger();

NullLogger& LoggerFactory::get_null_logger()
{
    return *null_logger_instance;
}
} // namespace Loggers