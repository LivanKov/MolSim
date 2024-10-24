#ifndef GETOPT_WRAPPER
#define GETOPT_WRAPPER

#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <string>
#include <unordered_map>
#include <string.h>

class GetoptWrapper {
public:
    //strdup helps negate compiler warnings
  GetoptWrapper(const std::string_view &optstring)
      : optstring_{optstring}, opt_map_{{'e', strdup("1000")},
                                        {'d', strdup("0.014")},
                                        {'i', strdup("../input/eingabe-sonne.txt")},
                                        {'o', strdup("../output")},
                                        {'s', strdup("1")}} {}
  // initialize the map with default values
  // default values can be changed later on
  std::unordered_map<char, char *> &parse(int argc, char **argv) {
    int opt;
    while ((opt = getopt(argc, argv, optstring_.c_str())) != -1) {
      switch (opt) {
      case 'e':
        opt_map_['e'] = optarg;
        break;
      case 'd':
        opt_map_['d'] = optarg;
        break;
      case 'i':
        opt_map_['i'] = optarg;
        break;
      case 'o':
        opt_map_['o'] = optarg;
        break;
      case 's':
        opt_map_['s'] = optarg;
        break;
      default:
        std::cerr << "Unknown option or missing argument.\n";
        break;
      }
    }
    return opt_map_;
  }

private:
  std::string optstring_;
  std::unordered_map<char, char*> opt_map_;
};

#endif
