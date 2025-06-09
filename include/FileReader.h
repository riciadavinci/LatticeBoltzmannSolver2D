#pragma once
#include <fstream>
#include <map>
#include <string>
#include <sstream>
#include <iostream>

class FileReader
{
private:
   std::map<std::string, std::string> data_;

public:
   FileReader(std::string filename);

   template<typename T>
   T getParam(const std::string& key);
};

#include "../src/FileReader.cpp"
