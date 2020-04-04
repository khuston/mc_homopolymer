#include <string>
#include <memory>
#include <vector>

#ifndef LOGGERS_H
#define LOGGERS_H

namespace Loggers
{
enum LogLevel
{
    Debug,
    Info,
    Warn,
    Error
};
enum LoggerType
{
    Stdout
};

class Logger
{
  protected:
    friend class LoggerFactory;
    LogLevel level;
    bool should_log (const LogLevel& message_level) const;

  public:
    Logger();
    virtual void log (const LogLevel& message_level, const std::string& message) const = 0;
};

class StdoutLogger : public Logger
{
  public:
    using Logger::Logger;
    void log (const LogLevel& message_level, const std::string& message) const override;
};

class CombinedLogger : public Logger
{
  private:
    std::unique_ptr<std::vector<Logger>> loggers;
  public:
    CombinedLogger();
    void log (const LogLevel& message_level, const std::string& message) const override;

};

class LoggerFactory
{
  public:
    static std::unique_ptr<Logger> make_stdout_logger(const LogLevel& level);
    static std::unique_ptr<Logger> make_combined_logger(std::initializer_list<Logger> loggers_to_combine);
};
} // namespace Loggers

#endif