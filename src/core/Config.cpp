#include "core/Config.h"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace jrc {

std::size_t dataTypeBytes(DataType type)
{
    switch (type) {
        case DataType::Int16:
        case DataType::Uint16: return 2;
        case DataType::Int32:
        case DataType::Uint32:
        case DataType::Single: return 4;
        case DataType::Double: return 8;
    }
    return 0;
}

DataType dataTypeFromString(const std::string& name)
{
    if (name == "int16") return DataType::Int16;
    if (name == "uint16") return DataType::Uint16;
    if (name == "int32") return DataType::Int32;
    if (name == "uint32") return DataType::Uint32;
    if (name == "single" || name == "float" || name == "float32") return DataType::Single;
    if (name == "double" || name == "float64") return DataType::Double;
    throw std::runtime_error("unsupported dataType: '" + name + "'");
}

const char* dataTypeName(DataType type)
{
    switch (type) {
        case DataType::Int16: return "int16";
        case DataType::Uint16: return "uint16";
        case DataType::Int32: return "int32";
        case DataType::Uint32: return "uint32";
        case DataType::Single: return "single";
        case DataType::Double: return "double";
    }
    return "?";
}

namespace {

std::string trim(const std::string& s)
{
    const auto begin = s.find_first_not_of(" \t\r\n");
    if (begin == std::string::npos) return {};
    const auto end = s.find_last_not_of(" \t\r\n");
    return s.substr(begin, end - begin + 1);
}

/// Removes `%` comments while respecting single-quoted strings, and drops the
/// MATLAB line-continuation token `...` so a wrapped array reads as one line.
std::string stripComments(const std::string& text)
{
    std::string out;
    out.reserve(text.size());

    bool inQuote = false;
    for (std::size_t i = 0; i < text.size(); ++i) {
        const char c = text[i];

        if (inQuote) {
            out.push_back(c);
            // '' inside a quoted string is an escaped quote, not a terminator.
            if (c == '\'') {
                if (i + 1 < text.size() && text[i + 1] == '\'') {
                    out.push_back(text[++i]);
                } else {
                    inQuote = false;
                }
            }
            continue;
        }

        if (c == '\'') {
            inQuote = true;
            out.push_back(c);
            continue;
        }

        if (c == '%') {  // comment runs to end of line
            while (i < text.size() && text[i] != '\n') ++i;
            out.push_back('\n');
            continue;
        }

        if (c == '.' && i + 2 < text.size() && text[i + 1] == '.' && text[i + 2] == '.') {
            i += 2;                                              // skip '...'
            while (i + 1 < text.size() && text[i + 1] != '\n') ++i;  // and the newline
            ++i;
            out.push_back(' ');
            continue;
        }

        out.push_back(c);
    }

    return out;
}

/// Splits `body` into statements at top-level `;` or newline.
std::vector<std::string> splitStatements(const std::string& body)
{
    std::vector<std::string> statements;
    std::string current;

    int depth = 0;
    bool inQuote = false;

    for (std::size_t i = 0; i < body.size(); ++i) {
        const char c = body[i];

        if (inQuote) {
            current.push_back(c);
            if (c == '\'') {
                if (i + 1 < body.size() && body[i + 1] == '\'') {
                    current.push_back(body[++i]);
                } else {
                    inQuote = false;
                }
            }
            continue;
        }

        switch (c) {
            case '\'': inQuote = true; current.push_back(c); break;
            case '[':
            case '{':
            case '(': ++depth; current.push_back(c); break;
            case ']':
            case '}':
            case ')': --depth; current.push_back(c); break;
            case ';':
            case '\n':
                if (depth <= 0) {
                    if (!trim(current).empty()) statements.push_back(trim(current));
                    current.clear();
                } else if (c != '\n') {
                    current.push_back(c);
                } else {
                    current.push_back(' ');
                }
                break;
            default: current.push_back(c);
        }
    }

    if (!trim(current).empty()) statements.push_back(trim(current));
    return statements;
}

/// Parses a right-hand side into numbers and/or strings.
ParamValue parseValue(const std::string& rhsIn)
{
    ParamValue value;
    value.raw = trim(rhsIn);

    std::string rhs = value.raw;
    if (rhs.empty()) return value;

    // Peel one layer of [] or {}; both hold comma/space separated elements.
    if ((rhs.front() == '[' && rhs.back() == ']') ||
        (rhs.front() == '{' && rhs.back() == '}')) {
        rhs = trim(rhs.substr(1, rhs.size() - 2));
    }
    if (rhs.empty()) return value;  // [] or {}

    // Tokenise on commas and whitespace, keeping quoted strings intact.
    std::vector<std::string> tokens;
    std::string current;
    bool inQuote = false;

    for (std::size_t i = 0; i < rhs.size(); ++i) {
        const char c = rhs[i];

        if (inQuote) {
            if (c == '\'') {
                if (i + 1 < rhs.size() && rhs[i + 1] == '\'') {
                    current.push_back('\'');
                    ++i;
                } else {
                    inQuote = false;
                    tokens.push_back(current);
                    current.clear();
                }
            } else {
                current.push_back(c);
            }
            continue;
        }

        if (c == '\'') {
            inQuote = true;
            current.clear();
            continue;
        }

        // ';' separates matrix rows; ',' and whitespace separate elements.
        if (c == ',' || c == ' ' || c == '\t' || c == ';') {
            if (!current.empty()) {
                tokens.push_back(current);
                current.clear();
            }
            if (c == ';') ++value.rows;
            continue;
        }

        current.push_back(c);
    }
    if (!current.empty()) tokens.push_back(current);

    // A token that parses cleanly as a number is numeric; otherwise it is text.
    const bool quoted = value.raw.find('\'') != std::string::npos;

    for (const auto& token : tokens) {
        if (quoted) {
            value.strings.push_back(token);
            continue;
        }
        try {
            std::size_t consumed = 0;
            const double number = std::stod(token, &consumed);
            if (consumed == token.size()) {
                value.numbers.push_back(number);
            } else {
                value.strings.push_back(token);
            }
        } catch (const std::exception&) {
            value.strings.push_back(token);
        }
    }

    return value;
}

bool isAbsolutePath(const std::string& path)
{
    if (path.empty()) return false;
    if (path[0] == '/' || path[0] == '\\') return true;
    return path.size() > 2 && path[1] == ':';  // C:\ or C:/
}

std::string directoryOf(const std::string& path)
{
    const auto pos = path.find_last_of("/\\");
    return pos == std::string::npos ? std::string(".") : path.substr(0, pos);
}

}  // namespace

void Config::load(const std::string& prmPath)
{
    std::ifstream stream(prmPath, std::ios::binary);
    if (!stream) throw std::runtime_error("cannot open parameter file: " + prmPath);

    std::ostringstream buffer;
    buffer << stream.rdbuf();

    path_ = prmPath;
    directory_ = directoryOf(prmPath);
    params_.clear();

    const std::string body = stripComments(buffer.str());

    for (const auto& statement : splitStatements(body)) {
        // Split at the first '=' that is not part of ==, <=, >= or ~=.
        std::size_t eq = std::string::npos;
        int depth = 0;
        bool inQuote = false;

        for (std::size_t i = 0; i < statement.size(); ++i) {
            const char c = statement[i];
            if (inQuote) {
                if (c == '\'') inQuote = false;
                continue;
            }
            if (c == '\'') { inQuote = true; continue; }
            if (c == '[' || c == '{' || c == '(') { ++depth; continue; }
            if (c == ']' || c == '}' || c == ')') { --depth; continue; }
            if (c == '=' && depth == 0) {
                const bool comparison =
                    (i + 1 < statement.size() && statement[i + 1] == '=') ||
                    (i > 0 && (statement[i - 1] == '=' || statement[i - 1] == '<' ||
                               statement[i - 1] == '>' || statement[i - 1] == '~' ||
                               statement[i - 1] == '!'));
                if (!comparison) { eq = i; break; }
            }
        }

        if (eq == std::string::npos) continue;  // not an assignment

        const std::string key = trim(statement.substr(0, eq));
        const std::string rhs = statement.substr(eq + 1);
        if (key.empty()) continue;

        params_[key] = parseValue(rhs);
    }

    resolveRecordingParams();
}

void Config::resolveRecordingParams()
{
    nChans_ = getInt("nChans", 0);
    sampleRate_ = getNumber("sampleRate", 0.0);
    headerOffset_ = static_cast<std::size_t>(getNumber("headerOffset", 0.0));
    bitScaling_ = getNumber("bitScaling", 1.0);

    const std::string typeName = getString("dataType", "int16");
    dataType_ = dataTypeFromString(typeName);

    rawRecordings_.clear();
    for (const auto& entry : getStrings("rawRecordings")) {
        if (entry.empty()) continue;
        rawRecordings_.push_back(isAbsolutePath(entry) ? entry
                                                       : directory_ + "/" + entry);
    }
}

bool Config::has(const std::string& key) const
{
    return params_.find(key) != params_.end();
}

const ParamValue* Config::find(const std::string& key) const
{
    const auto it = params_.find(key);
    return it == params_.end() ? nullptr : &it->second;
}

double Config::getNumber(const std::string& key, double fallback) const
{
    const auto* value = find(key);
    return (value && !value->numbers.empty()) ? value->numbers.front() : fallback;
}

int Config::getInt(const std::string& key, int fallback) const
{
    const auto* value = find(key);
    return (value && !value->numbers.empty())
               ? static_cast<int>(value->numbers.front())
               : fallback;
}

std::string Config::getString(const std::string& key, const std::string& fallback) const
{
    const auto* value = find(key);
    return (value && !value->strings.empty()) ? value->strings.front() : fallback;
}

std::vector<double> Config::getNumbers(const std::string& key) const
{
    const auto* value = find(key);
    return value ? value->numbers : std::vector<double>{};
}

std::vector<std::string> Config::getStrings(const std::string& key) const
{
    const auto* value = find(key);
    return value ? value->strings : std::vector<std::string>{};
}

std::vector<double> Config::getMatrix(const std::string& key, std::size_t columns) const
{
    const auto* value = find(key);
    if (!value || value->numbers.empty()) return {};

    if (columns == 0 || value->numbers.size() % columns != 0) {
        throw std::runtime_error("parameter '" + key + "' holds " +
                                 std::to_string(value->numbers.size()) +
                                 " values, which is not a multiple of " +
                                 std::to_string(columns));
    }

    return value->numbers;
}

}  // namespace jrc
