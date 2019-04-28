#ifndef SEQAN_HEADER_ARGS_PARSER_H
#define SEQAN_HEADER_ARGS_PARSER_H
#include "base.h"

seqan::ArgumentParser::ParseResult
parseCommandLine(Options & options, int argc, char const ** argv);

#endif
