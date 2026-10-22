//! \file regression_test.cpp
//! \brief Regression test for complete EPANET-CGGA simulations.
//!
//! The program runs one network from start to end with the adaptive solver
//! and compares the result with a reference file stored in the repository.
//! Three quantities are checked:
//!
//!   1. the head series at the monitored nodes (maximum absolute difference),
//!   2. the peak and minimum head at each monitored node,
//!   3. the mode schedule, i.e. the sequence of QS, RWC and WH intervals and
//!      the times at which the solver switches between them.
//!
//! The mode of each interval is identified from the time step size:
//! dt < 0.5 s is the water hammer mode, 0.5 s <= dt < 5 s is the rigid water
//! column mode, and dt >= 5 s is the quasi-steady mode.
//!
//! Usage
//!   regression_test --inp <file.inp> --ref <reference.csv> --out <result.csv>
//!                   --nodes <id1,id2,...> [--tol <m>] [--mode-tol <s>]
//!                   [--update]
//!
//! With --update the reference file is (re)written from the current run and
//! no comparison is made. Without it the program returns 0 if all checks
//! pass and 1 otherwise. The file given by --out is always written, so that
//! a failed comparison can be inspected.
//!
//! No external test framework is required.

#include "epanet3.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace {

// ---------------------------------------------------------------------------
// Data structures
// ---------------------------------------------------------------------------

struct Series
{
    std::vector<std::string>          nodes;
    std::vector<double>               time;   // s
    std::vector<std::vector<double>>  head;   // head[k][i]: node k, row i
};

struct ModeChange
{
    std::string mode;   // "QS", "RWC" or "WH"
    double      start;  // s
};

struct Extreme
{
    double value;
    double time;
};

struct Summary
{
    std::string              inpName;
    std::vector<std::string> nodes;
    std::vector<ModeChange>  modes;
    std::vector<Extreme>     peak;
    std::vector<Extreme>     minimum;
    long                     steps = 0;
};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

std::vector<std::string> split(const std::string& s, char sep)
{
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, sep))
    {
        // trim spaces and carriage returns
        size_t a = item.find_first_not_of(" \t\r");
        size_t b = item.find_last_not_of(" \t\r");
        if (a == std::string::npos) out.push_back("");
        else out.push_back(item.substr(a, b - a + 1));
    }
    return out;
}

std::string baseName(const std::string& path)
{
    size_t k = path.find_last_of("/\\");
    return (k == std::string::npos) ? path : path.substr(k + 1);
}

const char* modeOf(double dt)
{
    if (dt < 0.5) return "WH";
    if (dt < 5.0) return "RWC";
    return "QS";
}

// Mode schedule from the time series. The last interval is ignored because
// the final step is shortened to the end of the simulation.
std::vector<ModeChange> modeSchedule(const std::vector<double>& t)
{
    std::vector<ModeChange> out;
    for (size_t i = 0; i + 2 < t.size(); ++i)
    {
        const char* m = modeOf(t[i + 1] - t[i]);
        if (out.empty() || out.back().mode != m) out.push_back({ m, t[i] });
    }
    return out;
}

Summary summarize(const Series& s, const std::string& inpName)
{
    Summary sum;
    sum.inpName = inpName;
    sum.nodes   = s.nodes;
    sum.modes   = modeSchedule(s.time);
    sum.steps   = (long)s.time.size();
    for (size_t k = 0; k < s.nodes.size(); ++k)
    {
        const std::vector<double>& h = s.head[k];
        size_t iMax = std::max_element(h.begin(), h.end()) - h.begin();
        size_t iMin = std::min_element(h.begin(), h.end()) - h.begin();
        sum.peak.push_back({ h[iMax], s.time[iMax] });
        sum.minimum.push_back({ h[iMin], s.time[iMin] });
    }
    return sum;
}

// Linear interpolation of series column k at time t (time is increasing).
double interp(const Series& s, size_t k, double t)
{
    const std::vector<double>& x = s.time;
    const std::vector<double>& y = s.head[k];
    if (t <= x.front()) return y.front();
    if (t >= x.back())  return y.back();
    size_t j = std::upper_bound(x.begin(), x.end(), t) - x.begin();
    double x0 = x[j - 1], x1 = x[j];
    if (x1 - x0 <= 0.0) return y[j];
    double w = (t - x0) / (x1 - x0);
    return y[j - 1] + w * (y[j] - y[j - 1]);
}

// ---------------------------------------------------------------------------
// Simulation
// ---------------------------------------------------------------------------

int runSimulation(const std::string& inp, const std::string& rpt,
                  const std::vector<std::string>& nodes, Series& out)
{
    EN_Project p = EN_createProject();
    int err = 0;

    std::vector<int> index(nodes.size(), -1);
    out.nodes = nodes;
    out.head.assign(nodes.size(), std::vector<double>());

    do
    {
        if ((err = EN_openReportFile(rpt.c_str(), p))) break;
        if ((err = EN_loadProject(inp.c_str(), p)))    break;
        if ((err = EN_openOutputFile("", p)))          break;

        for (size_t k = 0; k < nodes.size(); ++k)
        {
            std::vector<char> name(nodes[k].begin(), nodes[k].end());
            name.push_back('\0');
            if (EN_getNodeIndex(name.data(), &index[k], p) != 0 || index[k] < 0)
            {
                std::fprintf(stderr, "Node %s not found in %s\n",
                             nodes[k].c_str(), inp.c_str());
                err = 1;
            }
        }
        if (err) break;

        if ((err = EN_initSolver(0, p))) break;

        double t = 0.0;
        double tstep = 0.0;
        do
        {
            if ((err = EN_runSolver(&t, p))) break;
            if ((err = EN_advanceSolver(&tstep, p))) break;

            out.time.push_back(t);
            for (size_t k = 0; k < nodes.size(); ++k)
            {
                double h = 0.0;
                EN_getNodeValue(index[k], EN_HEAD, &h, p);
                out.head[k].push_back(h);
            }
        } while (tstep > 0.0 && !err);
    } while (false);

    EN_deleteProject(p);
    return err;
}

// ---------------------------------------------------------------------------
// Reference file
// ---------------------------------------------------------------------------
//
// Format:
//   # EPANET-CGGA regression reference
//   # input: <file name>
//   # steps: <number of solver steps>
//   # modes: QS@0.000 WH@30.000 RWC@374.000 QS@380.000
//   # peak: J16=211.081@45.570 J24=...
//   # min: J16=176.900@165.430 J24=...
//   time,J16,J24
//   0.0000,195.871000,196.750000
//   ...
//
// Rows are stored at every step that falls on a multiple of 0.1 s, so that
// the file stays small. Peaks, minima and the mode schedule are taken from
// the full-resolution run.

bool onDecimatedGrid(double t)
{
    double r = t * 10.0;
    return std::fabs(r - std::round(r)) < 1e-6;
}

bool writeCsv(const std::string& path, const Series& s, const Summary& sum)
{
    std::ofstream f(path);
    if (!f) return false;
    f << "# EPANET-CGGA regression reference\n";
    f << "# input: " << sum.inpName << "\n";
    f << "# steps: " << sum.steps << "\n";
    f << "# modes:";
    char buf[128];
    for (const ModeChange& m : sum.modes)
    {
        std::snprintf(buf, sizeof(buf), " %s@%.3f", m.mode.c_str(), m.start);
        f << buf;
    }
    f << "\n# peak:";
    for (size_t k = 0; k < s.nodes.size(); ++k)
    {
        std::snprintf(buf, sizeof(buf), " %s=%.6f@%.3f", s.nodes[k].c_str(),
                      sum.peak[k].value, sum.peak[k].time);
        f << buf;
    }
    f << "\n# min:";
    for (size_t k = 0; k < s.nodes.size(); ++k)
    {
        std::snprintf(buf, sizeof(buf), " %s=%.6f@%.3f", s.nodes[k].c_str(),
                      sum.minimum[k].value, sum.minimum[k].time);
        f << buf;
    }
    f << "\ntime";
    for (const std::string& n : s.nodes) f << "," << n;
    f << "\n";
    for (size_t i = 0; i < s.time.size(); ++i)
    {
        if (!onDecimatedGrid(s.time[i])) continue;
        std::snprintf(buf, sizeof(buf), "%.4f", s.time[i]);
        f << buf;
        for (size_t k = 0; k < s.nodes.size(); ++k)
        {
            std::snprintf(buf, sizeof(buf), ",%.6f", s.head[k][i]);
            f << buf;
        }
        f << "\n";
    }
    return static_cast<bool>(f);
}

bool parseExtremes(const std::string& line, std::vector<Extreme>& out)
{
    // " J16=211.081000@45.570 J24=..."
    std::stringstream ss(line);
    std::string tok;
    while (ss >> tok)
    {
        size_t e = tok.find('='), a = tok.find('@');
        if (e == std::string::npos || a == std::string::npos) return false;
        out.push_back({ std::atof(tok.substr(e + 1, a - e - 1).c_str()),
                        std::atof(tok.substr(a + 1).c_str()) });
    }
    return true;
}

bool readCsv(const std::string& path, Series& s, Summary& sum)
{
    std::ifstream f(path);
    if (!f) return false;
    std::string line;
    bool header = false;
    while (std::getline(f, line))
    {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        if (line[0] == '#')
        {
            if (line.rfind("# input:", 0) == 0) sum.inpName = line.substr(9);
            else if (line.rfind("# steps:", 0) == 0) sum.steps = std::atol(line.substr(8).c_str());
            else if (line.rfind("# modes:", 0) == 0)
            {
                std::stringstream ss(line.substr(8));
                std::string tok;
                while (ss >> tok)
                {
                    size_t a = tok.find('@');
                    sum.modes.push_back({ tok.substr(0, a), std::atof(tok.substr(a + 1).c_str()) });
                }
            }
            else if (line.rfind("# peak:", 0) == 0) parseExtremes(line.substr(7), sum.peak);
            else if (line.rfind("# min:", 0) == 0)  parseExtremes(line.substr(6), sum.minimum);
            continue;
        }
        std::vector<std::string> cols = split(line, ',');
        if (!header)
        {
            header = true;
            s.nodes.assign(cols.begin() + 1, cols.end());
            s.head.assign(s.nodes.size(), std::vector<double>());
            sum.nodes = s.nodes;
            continue;
        }
        if (cols.size() != s.nodes.size() + 1) return false;
        s.time.push_back(std::atof(cols[0].c_str()));
        for (size_t k = 0; k < s.nodes.size(); ++k)
            s.head[k].push_back(std::atof(cols[k + 1].c_str()));
    }
    return header && !s.time.empty() && sum.peak.size() == s.nodes.size()
        && sum.minimum.size() == s.nodes.size() && !sum.modes.empty();
}

// ---------------------------------------------------------------------------
// Comparison
// ---------------------------------------------------------------------------

int failures = 0;

void check(bool ok, const std::string& what)
{
    std::printf("  %s  %s\n", ok ? "ok  " : "FAIL", what.c_str());
    std::fflush(stdout);
    if (!ok) ++failures;
}

std::string fmt(const char* f, double a, double b = 0.0, double c = 0.0, double d = 0.0)
{
    char buf[256];
    std::snprintf(buf, sizeof(buf), f, a, b, c, d);
    return buf;
}

void compare(const Series& ref, const Summary& refSum,
             const Series& cand, const Summary& candSum,
             double tol, double modeTol)
{
    // 1. mode schedule
    std::string refModes, candModes;
    for (const ModeChange& m : refSum.modes)  refModes  += " " + m.mode + fmt("@%.2f", m.start);
    for (const ModeChange& m : candSum.modes) candModes += " " + m.mode + fmt("@%.2f", m.start);
    std::printf("  reference modes:%s\n  current modes:  %s\n", refModes.c_str(), candModes.c_str());

    bool sameSequence = refSum.modes.size() == candSum.modes.size();
    double worstShift = 0.0;
    if (sameSequence)
    {
        for (size_t i = 0; i < refSum.modes.size(); ++i)
        {
            if (refSum.modes[i].mode != candSum.modes[i].mode) sameSequence = false;
            worstShift = std::max(worstShift,
                                  std::fabs(refSum.modes[i].start - candSum.modes[i].start));
        }
    }
    check(sameSequence, "mode sequence is unchanged");
    if (sameSequence)
        check(worstShift <= modeTol,
              fmt("mode transition times agree within %.2f s (largest shift %.2f s)", modeTol, worstShift));

    // 2. head series, peaks and minima
    for (size_t k = 0; k < ref.nodes.size(); ++k)
    {
        const std::string& node = ref.nodes[k];
        size_t kc = std::find(cand.nodes.begin(), cand.nodes.end(), node) - cand.nodes.begin();
        if (kc == cand.nodes.size()) { check(false, "node " + node + " present in current run"); continue; }

        double worst = 0.0, tWorst = 0.0;
        for (size_t i = 0; i < ref.time.size(); ++i)
        {
            double d = std::fabs(interp(cand, kc, ref.time[i]) - ref.head[k][i]);
            if (d > worst) { worst = d; tWorst = ref.time[i]; }
        }
        check(worst <= tol,
              "node " + node + fmt(": max |dH| = %.4f m at t = %.2f s (tolerance %.4f m)", worst, tWorst, tol));

        double dPeak = std::fabs(candSum.peak[kc].value - refSum.peak[k].value);
        check(dPeak <= tol,
              "node " + node + fmt(": peak head %.4f m, reference %.4f m", candSum.peak[kc].value, refSum.peak[k].value));

        double dMin = std::fabs(candSum.minimum[kc].value - refSum.minimum[k].value);
        check(dMin <= tol,
              "node " + node + fmt(": minimum head %.4f m, reference %.4f m", candSum.minimum[kc].value, refSum.minimum[k].value));
    }
}

void usage()
{
    std::printf("Usage: regression_test --inp <file.inp> --ref <reference.csv> --out <result.csv>\n"
                "                       --nodes <id1,id2,...> [--tol <m>] [--mode-tol <s>] [--update]\n");
}

} // namespace

int main(int argc, char* argv[])
{
    std::string inp, ref, out, nodeList;
    double tol = 0.01;      // m
    double modeTol = 1.0;   // s
    bool update = false;

    for (int i = 1; i < argc; ++i)
    {
        std::string a = argv[i];
        auto next = [&](void) -> std::string { return (i + 1 < argc) ? argv[++i] : ""; };
        if      (a == "--inp")      inp = next();
        else if (a == "--ref")      ref = next();
        else if (a == "--out")      out = next();
        else if (a == "--nodes")    nodeList = next();
        else if (a == "--tol")      tol = std::atof(next().c_str());
        else if (a == "--mode-tol") modeTol = std::atof(next().c_str());
        else if (a == "--update")   update = true;
        else { usage(); return 2; }
    }
    if (inp.empty() || ref.empty() || out.empty() || nodeList.empty()) { usage(); return 2; }

    std::vector<std::string> nodes = split(nodeList, ',');
    std::printf("Regression test: %s\n", baseName(inp).c_str());
    std::fflush(stdout);

    Series cand;
    int err = runSimulation(inp, out + ".rpt", nodes, cand);
    if (err || cand.time.size() < 3)
    {
        std::printf("  FAIL  simulation did not complete (error code %d)\n", err);
        return 1;
    }
    Summary candSum = summarize(cand, baseName(inp));
    std::printf("  simulation completed: %ld steps, end time %.2f s\n",
                candSum.steps, cand.time.back());
    if (!writeCsv(out, cand, candSum))
        std::printf("  warning: result file %s could not be written\n", out.c_str());

    if (update)
    {
        if (!writeCsv(ref, cand, candSum))
        {
            std::printf("  FAIL  reference file %s could not be written\n", ref.c_str());
            return 1;
        }
        std::printf("  reference written to %s\n", ref.c_str());
        return 0;
    }

    Series refSeries;
    Summary refSum;
    if (!readCsv(ref, refSeries, refSum))
    {
        std::printf("  FAIL  reference file %s is missing or unreadable.\n"
                    "        Generate it with --update (see tests/regression/README.md).\n",
                    ref.c_str());
        return 1;
    }

    compare(refSeries, refSum, cand, candSum, tol, modeTol);
    std::printf("%s\n", failures ? "REGRESSION TEST FAILED" : "Regression test passed");
    return failures ? 1 : 0;
}
