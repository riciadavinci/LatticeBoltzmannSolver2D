#include "FileReader.h"

FileReader::FileReader(std::string filename)
{
    std::ifstream istream(filename);
    std::string line;
    std::string key;
    std::string value;
    while (std::getline(istream, line))
    {
        std::stringstream ss(line);
        ss >> key >> value;
        data_.insert(std::make_pair(key, value));
    }
}

template<typename T>
T FileReader::getParam(const std::string& key)
{
    double param = std::stod(data_[key]);
    return static_cast<T>(param);
}

template<>
std::string FileReader::getParam<std::string>(const std::string& key)
{
    return data_[key];
}