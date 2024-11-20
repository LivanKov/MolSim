#ifndef FILEWRITER_H
#define FILEWRITER_H

#include <string>

namespace output{

class FileWriter {
public:
    virtual ~FileWriter() = default;

    virtual void write_file(const std::string& filepath, int iteration) = 0;

    //virtual void close() = 0;
};
}
#endif