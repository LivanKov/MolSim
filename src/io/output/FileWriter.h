#ifndef FILEWRITER_H
#define FILEWRITER_H

#include <string>

namespace output{

class FileWriter {
public:
    virtual ~FileWriter() = default;

    //virtual void write(const std::string& data) = 0;

    //virtual void close() = 0;
};
}
#endif