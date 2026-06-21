#pragma once
#include <string>
#include <vector>
#include <stdexcept>
#include <zlib.h>
#include <cstring>

namespace noLZSS {

/**
 * @brief RAII wrapper for writing gzipped binary data.
 * 
 * This class provides a convenient interface for writing binary data to gzipped files
 * using the zlib library. It automatically handles compression and resource cleanup.
 */
class GzipWriter {
private:
    gzFile file_;
    std::string path_;
    bool is_open_;

public:
    /**
     * @brief Opens a file for gzipped writing.
     * 
     * @param path Output file path (typically with .gz extension)
     * @param compression_level Compression level (0-9, default 6)
     * @throws std::runtime_error If file cannot be opened
     */
    explicit GzipWriter(const std::string& path, int compression_level = 6)
        : path_(path), is_open_(false) {
        
        // Open file with compression level
        std::string mode = "wb" + std::to_string(compression_level);
        file_ = gzopen(path.c_str(), mode.c_str());
        
        if (file_ == nullptr) {
            throw std::runtime_error("Cannot open gzip file for writing: " + path);
        }
        
        is_open_ = true;
    }

    /**
     * @brief Destructor - closes the file if still open.
     */
    ~GzipWriter() {
        close();
    }

    // Delete copy constructor and assignment operator
    GzipWriter(const GzipWriter&) = delete;
    GzipWriter& operator=(const GzipWriter&) = delete;

    /**
     * @brief Writes binary data to the gzipped file.
     * 
     * @param data Pointer to data to write
     * @param size Number of bytes to write
     * @throws std::runtime_error If write fails
     */
    void write(const void* data, size_t size) {
        if (!is_open_) {
            throw std::runtime_error("Attempt to write to closed gzip file: " + path_);
        }
        
        int written = gzwrite(file_, data, static_cast<unsigned int>(size));
        if (written != static_cast<int>(size)) {
            int errnum;
            const char* errmsg = gzerror(file_, &errnum);
            throw std::runtime_error("Failed to write to gzip file " + path_ + 
                                    ": " + std::string(errmsg));
        }
    }

    /**
     * @brief Flushes any buffered data to disk.
     * 
     * @throws std::runtime_error If flush fails
     */
    void flush() {
        if (is_open_) {
            if (gzflush(file_, Z_FINISH) != Z_OK) {
                int errnum;
                const char* errmsg = gzerror(file_, &errnum);
                throw std::runtime_error("Failed to flush gzip file " + path_ + 
                                        ": " + std::string(errmsg));
            }
        }
    }

    /**
     * @brief Closes the file.
     */
    void close() {
        if (is_open_) {
            gzclose(file_);
            is_open_ = false;
        }
    }

    /**
     * @brief Checks if the file is open.
     * 
     * @return true if the file is open, false otherwise
     */
    bool is_open() const {
        return is_open_;
    }
};

/**
 * @brief RAII wrapper for reading gzipped binary data.
 * 
 * This class provides a convenient interface for reading binary data from gzipped files.
 * It automatically detects whether a file is gzipped or not and handles decompression.
 */
class GzipReader {
private:
    gzFile file_;
    std::string path_;
    bool is_open_;

public:
    /**
     * @brief Opens a file for reading (auto-detects gzip compression).
     * 
     * @param path Input file path
     * @throws std::runtime_error If file cannot be opened
     */
    explicit GzipReader(const std::string& path)
        : path_(path), is_open_(false) {
        
        file_ = gzopen(path.c_str(), "rb");
        
        if (file_ == nullptr) {
            throw std::runtime_error("Cannot open file for reading: " + path);
        }
        
        is_open_ = true;
    }

    /**
     * @brief Destructor - closes the file if still open.
     */
    ~GzipReader() {
        close();
    }

    // Delete copy constructor and assignment operator
    GzipReader(const GzipReader&) = delete;
    GzipReader& operator=(const GzipReader&) = delete;

    /**
     * @brief Reads binary data from the file.
     * 
     * @param data Pointer to buffer to read into
     * @param size Number of bytes to read
     * @return Number of bytes actually read (may be less than size at EOF)
     * @throws std::runtime_error If read fails
     */
    size_t read(void* data, size_t size) {
        if (!is_open_) {
            throw std::runtime_error("Attempt to read from closed file: " + path_);
        }
        
        int bytes_read = gzread(file_, data, static_cast<unsigned int>(size));
        if (bytes_read < 0) {
            int errnum;
            const char* errmsg = gzerror(file_, &errnum);
            throw std::runtime_error("Failed to read from file " + path_ + 
                                    ": " + std::string(errmsg));
        }
        
        return static_cast<size_t>(bytes_read);
    }

    /**
     * @brief Seeks to a position in the uncompressed data stream.
     * 
     * @param offset Offset to seek to
     * @param whence SEEK_SET, SEEK_CUR, or SEEK_END
     * @return New position in the stream
     * @throws std::runtime_error If seek fails
     */
    z_off_t seek(z_off_t offset, int whence) {
        if (!is_open_) {
            throw std::runtime_error("Attempt to seek in closed file: " + path_);
        }
        
        z_off_t pos = gzseek(file_, offset, whence);
        if (pos < 0) {
            int errnum;
            const char* errmsg = gzerror(file_, &errnum);
            throw std::runtime_error("Failed to seek in file " + path_ + 
                                    ": " + std::string(errmsg));
        }
        
        return pos;
    }

    /**
     * @brief Gets the current position in the uncompressed data stream.
     * 
     * @return Current position
     */
    z_off_t tell() {
        if (!is_open_) {
            throw std::runtime_error("Attempt to tell on closed file: " + path_);
        }
        
        return gztell(file_);
    }

    /**
     * @brief Closes the file.
     */
    void close() {
        if (is_open_) {
            gzclose(file_);
            is_open_ = false;
        }
    }

    /**
     * @brief Checks if the file is open.
     * 
     * @return true if the file is open, false otherwise
     */
    bool is_open() const {
        return is_open_;
    }

    /**
     * @brief Checks if we've reached the end of file.
     * 
     * @return true if at EOF, false otherwise
     */
    bool eof() const {
        if (!is_open_) {
            return true;
        }
        return gzeof(file_);
    }
};

/**
 * @brief Checks if a file is gzip-compressed by reading its magic bytes.
 * 
 * @param path File path to check
 * @return true if the file is gzipped, false otherwise
 */
inline bool is_gzipped_file(const std::string& path) {
    // Open file in binary mode
    FILE* f = fopen(path.c_str(), "rb");
    if (f == nullptr) {
        return false;
    }
    
    // Read first two bytes (gzip magic: 0x1f 0x8b)
    unsigned char magic[2];
    size_t bytes_read = fread(magic, 1, 2, f);
    fclose(f);
    
    if (bytes_read != 2) {
        return false;
    }
    
    return (magic[0] == 0x1f && magic[1] == 0x8b);
}

} // namespace noLZSS
