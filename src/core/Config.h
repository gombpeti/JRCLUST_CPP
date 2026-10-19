#pragma once

#include <map>
#include <string>
#include <vector>

namespace jrc {

/// Sample format of the raw recording, from the `dataType` parameter.
enum class DataType { Int16, Uint16, Int32, Uint32, Single, Double };

/// Bytes occupied by one sample of `type`.
std::size_t dataTypeBytes(DataType type);

/// Parses a dataType string ('int16', 'single', ...). Throws if unrecognised.
DataType dataTypeFromString(const std::string& name);

const char* dataTypeName(DataType type);

/// One right-hand side from a .prm file.
///
/// Parameter files hold scalars, quoted strings, numeric arrays and cell arrays
/// of strings, so both representations are carried; the parser fills whichever
/// the syntax produced.
struct ParamValue {
    std::vector<double> numbers;
    std::vector<std::string> strings;
    std::string raw;  ///< text as written, for diagnostics

    /// Row count for matrix literals: `[a, b; c, d]` yields 2. Plain vectors
    /// and scalars yield 1. Values are stored row-major in `numbers`.
    std::size_t rows = 1;

    bool empty() const { return numbers.empty() && strings.empty(); }
};

/// A JRCLUST .prm parameter file.
///
/// Mirrors the subset of +jrclust/@Config needed to open a recording. Unknown
/// parameters are kept verbatim so nothing is silently dropped.
class Config {
public:
    /// Reads and parses `prmPath`. Throws std::runtime_error if unreadable.
    void load(const std::string& prmPath);

    bool has(const std::string& key) const;
    const ParamValue* find(const std::string& key) const;

    // -- typed accessors, each falling back to `fallback` when absent --------
    double getNumber(const std::string& key, double fallback) const;
    int getInt(const std::string& key, int fallback) const;
    std::string getString(const std::string& key, const std::string& fallback) const;
    std::vector<double> getNumbers(const std::string& key) const;
    std::vector<std::string> getStrings(const std::string& key) const;

    /// Reads an N x `columns` matrix, returned row-major. Throws if the value
    /// is present but its element count is not a multiple of `columns`.
    std::vector<double> getMatrix(const std::string& key, std::size_t columns) const;

    // -- recording geometry, resolved during load() -------------------------
    int nChans() const { return nChans_; }
    double sampleRate() const { return sampleRate_; }
    DataType dataType() const { return dataType_; }
    std::size_t headerOffset() const { return headerOffset_; }
    double bitScaling() const { return bitScaling_; }

    /// Raw recording paths from `rawRecordings`, resolved relative to the .prm.
    const std::vector<std::string>& rawRecordings() const { return rawRecordings_; }

    /// Replaces the recording list, for when a .prm carries stale paths
    /// (a drive letter that has since changed, say) and the data lives
    /// elsewhere. Leaves the file on disk untouched.
    void overrideRecordings(const std::vector<std::string>& paths)
    {
        rawRecordings_ = paths;
    }

    const std::string& path() const { return path_; }
    const std::map<std::string, ParamValue>& all() const { return params_; }

private:
    void resolveRecordingParams();

    std::string path_;
    std::string directory_;
    std::map<std::string, ParamValue> params_;

    int nChans_ = 0;
    double sampleRate_ = 0.0;
    DataType dataType_ = DataType::Int16;
    std::size_t headerOffset_ = 0;
    double bitScaling_ = 1.0;
    std::vector<std::string> rawRecordings_;
};

}  // namespace jrc
