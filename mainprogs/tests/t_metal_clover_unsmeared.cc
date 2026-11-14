// Stream a source-t=0 to unsmeared sink-t=0 propagator using Eloy's callback.
#include "actions/ferm/fermacts/eoprec_clover_fermact_w.h"
#include "actions/ferm/fermstates/ferm_createstate_aggregate_w.h"
#include "actions/ferm/invert/syssolver_linop_aggregate.h"
#include "actions/ferm/linop/unprec_clover_linop_w.h"
#include "util/ferm/mgproton.h"
#include "util/gauge/reunit.h"
#include <algorithm>
#include <bit>
#include <cerrno>
#include <cmath>
#include <csignal>
#include <cstdio>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <iomanip>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <sys/stat.h>
#include <sys/statvfs.h>
#include <unistd.h>
#include <vector>

#ifndef BUILD_SB
#error This program requires SuperbBLAS source/IO helpers and MetalClover
#endif
using namespace Chroma;
namespace {
  using T = LatticeFermion;
  using Q = multi1d<LatticeColorMatrix>;
  struct UnsmearedComplex {
    double re, im;
  };
  static_assert(sizeof(UnsmearedComplex) == 16);
  static_assert(std::endian::native == std::endian::little);
  struct AlreadyComplete {};
  volatile sig_atomic_t interrupted = 0;
  void on_signal(int) { interrupted = 1; }
  constexpr uint64_t fnv_offset = 14695981039346656037ull;
  uint64_t hash_bytes(const void *ptr, size_t size, uint64_t h = fnv_offset) {
    const auto *p = static_cast<const unsigned char *>(ptr);
    for (size_t i = 0; i < size; ++i)
      h = (h ^ p[i]) * 1099511628211ull;
    return h;
  }
  std::string hex64(uint64_t value) {
    std::ostringstream s;
    s << std::hex << std::setfill('0') << std::setw(16) << value;
    return s.str();
  }
  std::string read_text(const std::string &path) {
    std::ifstream f(path, std::ios::binary);
    if (!f)
      throw std::runtime_error("cannot read " + path);
    return {std::istreambuf_iterator<char>(f),
            std::istreambuf_iterator<char>()};
  }
  // QDP's string writer emits XML literally, so embedded documents must not
  // carry their own XML declaration inside the enclosing metadata document.
  std::string xml_fragment(std::string text) {
    const auto first = text.find_first_not_of(" \t\r\n");
    if (first != std::string::npos && text.compare(first, 5, "<?xml") == 0) {
      const auto end = text.find("?>", first);
      if (end == std::string::npos)
        throw std::runtime_error("unterminated XML declaration");
      text.erase(first, end + 2 - first);
    }
    return text;
  }
  std::string file_hash(const std::string &path) {
    std::ifstream f(path, std::ios::binary);
    if (!f)
      throw std::runtime_error("cannot hash " + path);
    std::vector<char> b(1 << 20);
    uint64_t h = fnv_offset;
    while (f) {
      f.read(b.data(), b.size());
      h = hash_bytes(b.data(), size_t(f.gcount()), h);
    }
    if (!f.eof())
      throw std::runtime_error("input hash read failed");
    return hex64(h);
  }
  bool exists(const std::string &name) {
    struct stat s{};
    return stat(name.c_str(), &s) == 0;
  }
  void io_error(const std::string &what) {
    throw std::runtime_error(what + ": " + std::strerror(errno));
  }
  void sync_fd(int fd) {
    if (fsync(fd))
      io_error("fsync");
  }
  void transfer(int fd, void *ptr, size_t count, uint64_t offset,
                bool writing) {
    auto *p = static_cast<char *>(ptr);
    while (count) {
      ssize_t n = writing ? pwrite(fd, p, count, off_t(offset))
                          : pread(fd, p, count, off_t(offset));
      if (n < 0 && errno == EINTR)
        continue;
      if (n <= 0)
        io_error(writing ? "payload write" : "payload read");
      p += n;
      count -= size_t(n);
      offset += uint64_t(n);
    }
  }
  std::string npy_header(int vectors, const multi1d<int> &dims) {
    std::ostringstream s;
    s << "{'descr': '<c16', 'fortran_order': False, 'shape': (" << vectors
      << ", 4, " << dims[2] << ", " << dims[1] << ", " << dims[0]
      << ", 3, 4), }";
    std::string h = s.str();
    h.append((64 - ((10 + h.size() + 1) % 64)) % 64, ' ');
    h += '\n';
    if (h.size() > 65535)
      throw std::runtime_error("NPY header too large");
    std::string out("\x93NUMPY\x01\x00", 8);
    out += char(h.size() & 255);
    out += char((h.size() >> 8) & 255);
    out += h;
    return out;
  }
  struct Payload {
    std::string prefix, signature, header;
    int fd = -1;
    FILE *csv = nullptr;
    uint64_t rows, total, completed = 0;
    double max_relative = 0, qprop_seconds = 0, verify_seconds = 0,
           write_seconds = 0;
    std::vector<UnsmearedComplex> readback;
    Payload(std::string p, std::string metadata, int vectors,
            const multi1d<int> &dims, bool resume)
        : prefix(std::move(p)), header(npy_header(vectors, dims)),
          rows(uint64_t(dims[0]) * dims[1] * dims[2] * Nc * Ns),
          total(uint64_t(vectors) * Ns), readback(rows) {
      signature = hex64(hash_bytes(metadata.data(), metadata.size()));
      uint64_t bytes = header.size() + rows * total * sizeof(UnsmearedComplex);
      if (resume) {
        if (read_text(prefix + ".meta.xml") != metadata)
          throw std::runtime_error("resume metadata/input mismatch");
        std::ifstream state(prefix + ".state");
        std::string magic, saved;
        uint64_t csv_bytes;
        if (!(state >> magic >> saved >> completed >> csv_bytes >>
              max_relative >> qprop_seconds >> verify_seconds >>
              write_seconds) ||
            magic != "METAL_UNSMEARED_V1" || saved != signature ||
            completed > total)
          throw std::runtime_error("invalid resume checkpoint");
        fd = open((prefix + ".npy").c_str(), O_RDWR);
        if (fd < 0)
          io_error("open resume payload");
        struct stat s{};
        if (fstat(fd, &s) || uint64_t(s.st_size) != bytes)
          throw std::runtime_error("resume payload size mismatch");
        std::string stored(header.size(), '\0');
        transfer(fd, stored.data(), stored.size(), 0, false);
        if (stored != header)
          throw std::runtime_error("resume NPY header mismatch");
        csv = fopen((prefix + ".columns.csv").c_str(), "r+b");
        if (!csv)
          io_error("open resume journal");
        if (ftruncate(fileno(csv), off_t(csv_bytes)))
          io_error("truncate uncommitted journal tail");
        if (fseeko(csv, 0, SEEK_END))
          io_error("seek journal");
        verify_journal(vectors);
      } else {
        for (const auto &suffix :
             {".npy", ".meta.xml", ".columns.csv", ".state"})
          if (exists(prefix + suffix))
            throw std::runtime_error(
                "output already exists; use --resume or a fresh prefix");
        fd = open((prefix + ".npy").c_str(), O_CREAT | O_EXCL | O_RDWR, 0644);
        if (fd < 0)
          io_error("create payload");
        fstore_t allocation{};
        allocation.fst_flags = F_ALLOCATEALL;
        allocation.fst_posmode = F_PEOFPOSMODE;
        allocation.fst_length = off_t(bytes);
        if (fcntl(fd, F_PREALLOCATE, &allocation))
          io_error("reserve output disk space");
        if (ftruncate(fd, off_t(bytes)))
          io_error("size payload");
        transfer(fd, header.data(), header.size(), 0, true);
        std::ofstream meta(prefix + ".meta.xml", std::ios::binary);
        meta << metadata;
        meta.close();
        if (!meta)
          throw std::runtime_error("write metadata failed");
        csv = fopen((prefix + ".columns.csv").c_str(), "w+b");
        if (!csv)
          io_error("create journal");
        fprintf(csv, "ordered_rhs,vector,source_spin,column,iterations,full_"
                     "relative,qprop_relative,source_norm,fnv1a64\n");
        checkpoint();
      }
    }
    ~Payload() {
      if (csv)
        fclose(csv);
      if (fd >= 0)
        close(fd);
    }
    void verify_journal(int vectors) {
      std::ifstream f(prefix + ".columns.csv");
      std::string line;
      std::getline(f, line);
      for (uint64_t i = 0; i < completed; ++i) {
        if (!std::getline(f, line))
          throw std::runtime_error("short checkpoint journal");
        std::vector<std::string> fields;
        std::istringstream row(line);
        std::string field;
        while (std::getline(row, field, ','))
          fields.push_back(field);
        if (fields.size() != 9 || std::stoull(fields[0]) != i)
          throw std::runtime_error("checkpoint journal order mismatch");
        uint64_t column = std::stoull(fields[3]);
        if (column != (i % vectors) * Ns + i / vectors)
          throw std::runtime_error("checkpoint column mismatch");
        if (!(std::stod(fields[5]) <= 1e-12 && std::stod(fields[6]) <= 1e-12))
          throw std::runtime_error("checkpoint tolerance mismatch");
        transfer(fd, readback.data(), rows * sizeof(UnsmearedComplex),
                 header.size() + column * rows * sizeof(UnsmearedComplex),
                 false);
        if (hex64(hash_bytes(readback.data(),
                             rows * sizeof(UnsmearedComplex))) != fields[8])
          throw std::runtime_error("checkpoint payload checksum mismatch");
      }
      std::cout << "RESUME_VERIFIED completed_rhs=" << completed << std::endl;
    }
    void save(int global_vector, int local_vector, int spin,
              const std::vector<UnsmearedComplex> &column, int iterations,
              double relative, double qrelative, double source_norm) {
      if (column.size() != rows)
        throw std::runtime_error("payload row count mismatch");
      uint64_t index = uint64_t(local_vector) * Ns + spin;
      uint64_t offset = header.size() + index * rows * sizeof(UnsmearedComplex);
      auto h =
          hash_bytes(column.data(), column.size() * sizeof(UnsmearedComplex));
      transfer(fd, const_cast<UnsmearedComplex *>(column.data()),
               column.size() * sizeof(UnsmearedComplex), offset, true);
      transfer(fd, readback.data(), column.size() * sizeof(UnsmearedComplex),
               offset, false);
      if (std::memcmp(readback.data(), column.data(),
                      column.size() * sizeof(UnsmearedComplex)))
        throw std::runtime_error("payload readback mismatch");
      fprintf(csv, "%llu,%d,%d,%llu,%d,%.17g,%.17g,%.17g,%s\n",
              (unsigned long long)completed, global_vector, spin,
              (unsigned long long)index, iterations, relative, qrelative,
              source_norm, hex64(h).c_str());
      if (ferror(csv))
        io_error("write column journal");
      max_relative = std::max(max_relative, relative);
      ++completed;
    }
    void checkpoint() {
      sync_fd(fd);
      if (fflush(csv))
        io_error("flush journal");
      sync_fd(fileno(csv));
      auto csv_bytes = ftello(csv);
      if (csv_bytes < 0)
        io_error("journal size");
      std::string temp = prefix + ".state.tmp";
      {
        std::ofstream f(temp, std::ios::binary | std::ios::trunc);
        f << std::setprecision(17) << "METAL_UNSMEARED_V1\n"
          << signature << '\n'
          << completed << '\n'
          << csv_bytes << '\n'
          << max_relative << '\n'
          << qprop_seconds << '\n'
          << verify_seconds << '\n'
          << write_seconds << '\n';
        f.close();
        if (!f)
          throw std::runtime_error("write checkpoint failed");
      }
      int state_fd = open(temp.c_str(), O_RDONLY);
      if (state_fd < 0)
        io_error("open checkpoint");
      sync_fd(state_fd);
      close(state_fd);
      if (rename(temp.c_str(), (prefix + ".state").c_str()))
        io_error("publish checkpoint");
    }
  };
  struct AuditedQprop : SystemSolver<T> {
    Handle<SystemSolver<T>> inner;
    mutable std::vector<SystemSolverResults_t> last;
    mutable std::vector<double> norms;
    mutable double last_seconds = 0;
    explicit AuditedQprop(Handle<SystemSolver<T>> p) : inner(p) {}
    const Subset &subset() const override { return inner->subset(); }
    SystemSolverResults_t operator()(T &x, const T &b) const override {
      return (*inner)(x, b);
    }
    std::vector<SystemSolverResults_t>
    operator()(const std::vector<std::shared_ptr<T>> &x,
               const std::vector<std::shared_ptr<const T>> &b) const override {
      StopWatch watch;
      watch.start();
      last = (*inner)(x, b);
      watch.stop();
      last_seconds = watch.getTimeInSeconds();
      if (last.size() != b.size())
        throw std::runtime_error("qprop result count mismatch");
      norms.resize(b.size());
      for (size_t i = 0; i < b.size(); ++i) {
        norms[i] = std::sqrt(toDouble(norm2(*b[i])));
        if (!(norms[i] > 0 && std::isfinite(norms[i])) ||
            !(toDouble(last[i].resid) / norms[i] <= 1e-12))
          throw std::runtime_error("original qprop full FP64 residual failed");
      }
      return last;
    }
  };
  void inverter(XMLBufferWriter &xml, const CloverFermActParams &p,
                const std::string &method) {
    push(xml, "InvertParam");
    write(xml, "invType", method);
    write(xml, "RsdTarget", Real(1e-12));
    write(xml, "MaxIter", 65536);
    write(xml, "RHSBatch", 16);
    multi1d<int> blocks(4);
    blocks[0] = blocks[1] = blocks[2] = 3;
    blocks[3] = 8;
    write(xml, "Blocks", blocks);
    write(xml, "SecondFilter", 16);
    write(xml, "BottomIterations", 4);
    push(xml, "CloverParams");
    write(xml, "Mass", p.Mass);
    write(xml, "u0", p.u0);
    write(xml, "AnisoParam", p.anisoParam);
    write(xml, "clovCoeffR", p.clovCoeffR);
    write(xml, "clovCoeffT", p.clovCoeffT);
    pop(xml);
    pop(xml);
  }
} // namespace
int main(int argc, char **argv) {
  if (argc < 9) {
    std::cerr << "Usage: t_metal_clover_unsmeared GAUGE ACTION_XML EIG_MOD "
                 "PREFIX FIRST_VECTOR NUM_VECTORS MAX_RHS METHOD [--resume] "
                 "[--stop-after-batches N]\n";
    return 2;
  }
  std::string gauge = argv[1], action_path = argv[2], eig_path = argv[3],
              prefix = argv[4], method = argv[8];
  int first = std::stoi(argv[5]), count = std::stoi(argv[6]),
      cap = std::stoi(argv[7]), stop_after = 0;
  bool resume = false;
  for (int i = 9; i < argc; ++i) {
    if (std::string(argv[i]) == "--resume")
      resume = true;
    else if (std::string(argv[i]) == "--stop-after-batches" && i + 1 < argc)
      stop_after = std::stoi(argv[++i]);
    else
      throw std::runtime_error("unknown option");
  }
  if (first < 0 || count < 1 || first + count > 5184 || cap < 1 || cap > 64 ||
      stop_after < 0 ||
      (method != "MetalMGClover" && method != "MetalBiCGClover"))
    throw std::runtime_error("invalid production arguments");
  int qargc = 1;
  QDP_initialize(&qargc, &argv);
  int status = 0;
  try {
    static_assert(BASE_PRECISION == 64 && Nc == 3 && Ns == 4,
                  "requires FP64 three-colour Wilson Chroma");
    std::signal(SIGTERM, on_signal);
    std::signal(SIGINT, on_signal);
    XMLReader metadata(action_path);
    multi1d<int> dims;
    read(metadata, "/DBMetaData/lattSize", dims);
    if (dims.size() != 4 || dims[0] != 12 || dims[1] != 12 || dims[2] != 12 ||
        dims[3] != 96)
      throw std::runtime_error("this production harness targets 12^3x96");
    Layout::setLattSize(dims);
    Layout::create();
    if (Layout::numNodes() != 1)
      throw std::runtime_error("one MPI rank required");
    std::cout << std::setprecision(17);
    const std::string action_node =
        "/DBMetaData/Params/Propagator/FermionAction";
    CloverFermActParams params(metadata, action_node);
    XMLReader prop_metadata(metadata, "/DBMetaData/Params/Propagator");
    auto action_group = readXMLGroup(prop_metadata, "FermionAction", "FermAct");
    XMLBufferWriter inv;
    inverter(inv, params, method);
    GroupXML_t inv_group;
    inv_group.id = method;
    inv_group.xml = inv.str();
    inv_group.path = "/InvertParam";
    Q u(4);
    XMLReader file_xml, record_xml;
    QDPFileReader reader(file_xml, gauge, QDPIO_SERIAL);
    read(reader, record_xml, u);
    if (reader.bad())
      throw std::runtime_error("gauge LIME/checksum read failed");
    reader.close();
    for (int mu = 0; mu < 4; ++mu)
      reunit(u[mu]);
    auto eigen_store = SB::openColorvecStorage({eig_path});
    std::string eigen_metadata;
    if (eigen_store.mod)
      eigen_store.mod->getUserdata(eigen_metadata);
    else if (eigen_store.mod_d)
      eigen_store.mod_d->getUserdata(eigen_metadata);
    else
      eigen_metadata = eigen_store.s3t.metadata;
    std::istringstream eigen_stream(eigen_metadata);
    XMLReader eigen_xml(eigen_stream);
    int available = 0, time = 0;
    read(eigen_xml, "/MODMetaData/num_vecs", available);
    read(eigen_xml, "/MODMetaData/t_slice", time);
    if (available != 5184 || time != 0)
      throw std::runtime_error(
          "expected all 5184 supplied eigenvectors at t=0");
    int precision = colorvecMODPrecision(eigen_metadata);
    XMLBufferWriter output_meta;
    push(output_meta, "MetalCloverUnsmeared");
    write(output_meta, "format_version", 1);
    write(output_meta, "lattSize", dims);
    write(output_meta, "source_time", 0);
    write(output_meta, "sink_time", 0);
    write(output_meta, "first_vector", first);
    write(output_meta, "num_vectors", count);
    write(output_meta, "source_spins", Ns);
    write(output_meta, "method", method);
    write(output_meta, "tolerance", Double(1e-12));
    write(output_meta, "dtype", std::string("little-endian complex128"));
    write(output_meta, "axes_C_order",
          std::string("source_vector,source_spin,z,y,x,colour,sink_spin"));
    write(
        output_meta, "source_order",
        std::string("preserved MOD file order; no sorting or phase rotation"));
    write(output_meta, "input_eigen_precision_bits", precision);
    write(output_meta, "GaugeFile", gauge);
    write(output_meta, "GaugeFNV1a64", file_hash(gauge));
    write(output_meta, "ActionFile", action_path);
    write(output_meta, "ActionFNV1a64", file_hash(action_path));
    write(output_meta, "ActionXML", xml_fragment(read_text(action_path)));
    write(output_meta, "EigenvectorFile", eig_path);
    write(output_meta, "EigenvectorFNV1a64", file_hash(eig_path));
    write(output_meta, "EigenvectorXML", xml_fragment(eigen_metadata));
    write(output_meta, "InvertParamXML", xml_fragment(inv.str()));
    pop(output_meta);
    Payload output(prefix, output_meta.str(), count, dims, resume);
    if (output.completed == output.total) {
      std::cout << "PASS unsmeared output already complete rhs=" << output.total
                << std::endl;
      SB::closeColorvecStorage(eigen_store);
      throw AlreadyComplete{};
    }
    auto all_eigenvectors =
        SB::getColorvecs<SB::ComplexD>(eigen_store, u, 3, 0, 1, first + count,
                                       "cxyzXnt", SB::Coor<3>{{}}, SB::OnHost);
    auto eigenvectors =
        all_eigenvectors.kvslice_from_size({{'n', first}}, {{'n', count}});
    SB::closeColorvecStorage(eigen_store);
    CreateFermStateEnv::registerAll();
    EvenOddPrecCloverFermActEnv::registerAll();
    LinOpSysSolverEnv::registerAll();
    StopWatch setup;
    setup.start();
    SB::ChimeraSolver solver(action_group, inv_group, u);
    setup.stop();
    UnprecCloverLinOp full(solver.state, params);
    auto *audit = new AuditedQprop(solver.PP);
    solver.PP = Handle<SystemSolver<T>>(audit);
    std::cout << "PRODUCTION_SETUP seconds=" << setup.getTimeInSeconds()
              << " method=" << method
              << " hierarchy_instances=" << (method == "MetalMGClover" ? 1 : 0)
              << " vectors=" << count << " physical_rhs=" << output.total
              << " source_time=0 sink_time=0" << std::endl;
    std::vector<int> sites(size_t(dims[0]) * dims[1] * dims[2]);
    multi1d<int> coord(4);
    coord[3] = 0;
    for (size_t i = 0; i < sites.size(); ++i) {
      size_t q = i;
      for (int d = 0; d < 3; ++d) {
        coord[d] = q % dims[d];
        q /= dims[d];
      }
      sites[i] = Layout::linearSiteIndex(coord);
    }
    T x = zero, b = zero, r = zero;
    std::vector<UnsmearedComplex> column(sites.size() * Nc * Ns);
    int batches = 0;
    StopWatch overall;
    overall.start();
    while (output.completed < output.total) {
      int spin = int(output.completed / count),
          begin = int(output.completed % count),
          width = std::min(cap, count - begin);
      auto eigen_tile =
          eigenvectors.kvslice_from_size({{'n', begin}}, {{'n', width}});
      SB::Tensor<2, SB::Complex> spins("Ss", {Ns, 1}, SB::OnHost,
                                       SB::OnEveryoneReplicated);
      spins.set_zero();
      spins.set({spin, 0}, 1.);
      uint64_t prior = output.completed;
      const SB::doInversionOnTSliceFn<SB::Complex> callback =
          [&](SB::Tensor<Nd + 5, SB::Complex> field, int source_spin, int n) {
            if (source_spin != 0 || n < 0 || n >= width ||
                uint64_t(spin * count + begin + n) != output.completed ||
                audit->last.size() != size_t(width))
              throw std::runtime_error("callback ordering mismatch");
            StopWatch verify;
            verify.start();
            x = zero;
            b = zero;
            field.rename_dims({{'s', '!'}, {'S', 's'}})
                .copyTo(SB::asTensorView(x).toComplex());
            eigen_tile.kvslice_from_size({{'n', n}}, {{'n', 1}})
                .rename_dims({{'n', 's'}})
                .copyTo(SB::asTensorView(b).toComplex().kvslice_from_size(
                    {{'s', spin}, {'t', 0}}, {{'s', 1}, {'t', 1}}));
            full(r, x, PLUS);
            r -= b;
            double norm = std::sqrt(toDouble(norm2(b))),
                   relative = std::sqrt(toDouble(norm2(r))) / norm;
            double qrelative = toDouble(audit->last[n].resid) / audit->norms[n];
            if (!(norm > 0 && std::isfinite(norm) && relative <= 1e-12 &&
                  qrelative <= 1e-12))
              throw std::runtime_error(
                  "independent full original Chroma FP64 residual failed");
            verify.stop();
            output.verify_seconds += verify.getTimeInSeconds();
            StopWatch io;
            io.start();
            for (size_t i = 0; i < sites.size(); ++i)
              for (int c = 0; c < Nc; ++c)
                for (int s = 0; s < Ns; ++s) {
                  auto z = x.elem(sites[i]).elem(s).elem(c);
                  column[(i * Nc + c) * Ns + s] = {double(z.real()),
                                                   double(z.imag())};
                }
            output.save(first + begin + n, begin + n, spin, column,
                        audit->last[n].n_count, relative, qrelative, norm);
            io.stop();
            output.write_seconds += io.getTimeInSeconds();
          };
      SB::doInversion(solver, eigen_tile, 0, spins, width, callback);
      if (output.completed != prior + width)
        throw std::runtime_error("incomplete callback batch");
      output.qprop_seconds += audit->last_seconds;
      output.checkpoint();
      ++batches;
      std::cout << "PRODUCTION_PROGRESS completed_rhs=" << output.completed
                << " total_rhs=" << output.total << " source_spin=" << spin
                << " next_vector=" << begin + width << " caller_rhs=" << width
                << " max_full_relative=" << output.max_relative
                << " qprop_seconds=" << output.qprop_seconds
                << " verify_seconds=" << output.verify_seconds
                << " write_seconds=" << output.write_seconds << std::endl;
      if (interrupted || (stop_after && batches >= stop_after)) {
        if (output.completed < output.total)
          status = 75;
        break;
      }
    }
    overall.stop();
    std::cout << (status == 75 ? "PAUSED" : "PASS")
              << " unsmeared_t0 completed_rhs=" << output.completed
              << " total_rhs=" << output.total
              << " session_seconds=" << overall.getTimeInSeconds()
              << " qprop_seconds=" << output.qprop_seconds
              << " independent_FP64_seconds=" << output.verify_seconds
              << " write_seconds=" << output.write_seconds
              << " maximum_full_relative=" << output.max_relative << std::endl;
  } catch (const AlreadyComplete &) {
    status = 0;
  } catch (const std::exception &e) {
    std::cerr << "ERROR " << e.what() << std::endl;
    status = 1;
  } catch (const std::string &e) {
    std::cerr << "ERROR " << e << std::endl;
    status = 1;
  }
  QDP_finalize();
  return status;
}
