//! Read names: spike's reads take the input's own name shape, marked `SPIKE_`.
//!
//! Tools read a flowcell position out of a read's name. Picard MarkDuplicates
//! (raredisease's duplicate marker) splits it on `:` and, when there are 5 or
//! 7 fields, takes the last three as tile, x and y to find optical duplicates;
//! any other name gets no position, and Picard warns once
//! (docs/superpowers/plans/2026-10-04-read-names.md). Real names differ by
//! machine and dataset, so spike learns the shape from the input's own names
//! and writes its reads in it. The machine field becomes `SPIKE_` and the
//! read's internal name, so spike's reads can still be told apart.

/// What every one of spike's read names starts with.
pub const MARK: &str = "SPIKE_";

/// The shape of the input's read names.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub enum NameShape {
    /// Names Picard reads no position from; spike's reads get none either.
    #[default]
    Other,
    /// `MACHINE:<middle>:TILE:X:Y`, the middle being `RUN:FLOWCELL:LANE`
    /// (7 fields) or `LANE` (5 fields): the input's most common middle, the
    /// tiles seen with it (never empty), and the range of x and of y.
    Illumina {
        middle: String,
        tiles: Vec<u32>,
        x: (u32, u32),
        y: (u32, u32),
    },
}

/// How many sampled names fell in each class.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct NameClasses {
    pub seven: usize,
    pub five: usize,
    pub other: usize,
}

/// A name Picard reads a position from.
#[derive(Debug, PartialEq, Eq)]
struct Positioned<'a> {
    fields: usize,
    middle: &'a [u8],
    tile: u32,
    x: u32,
    y: u32,
}

/// A field of 1-9 ASCII digits, as a number; 9 digits always fit Picard's int.
fn number(field: &[u8]) -> Option<u32> {
    if !(1..=9).contains(&field.len()) || !field.iter().all(u8::is_ascii_digit) {
        return None;
    }
    std::str::from_utf8(field).ok()?.parse().ok()
}

/// The position Picard reads from `name`: 5 or 7 `:` fields, the last three
/// numbers. `None` for any other name.
fn positioned(name: &[u8]) -> Option<Positioned<'_>> {
    let fields: Vec<&[u8]> = name.split(|&b| b == b':').collect();
    let n = fields.len();
    if n != 5 && n != 7 {
        return None;
    }
    let (tile, x, y) = (number(fields[n - 3])?, number(fields[n - 2])?, number(fields[n - 1])?);
    let tail = fields[n - 3].len() + fields[n - 2].len() + fields[n - 1].len() + 3;
    Some(Positioned { fields: n, middle: &name[fields[0].len() + 1..name.len() - tail], tile, x, y })
}

/// One step of splitmix64: a well-mixed 64-bit value from `state`.
fn splitmix64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// FNV-1a, 64-bit.
fn fnv1a64(bytes: &[u8]) -> u64 {
    bytes.iter().fold(0xCBF2_9CE4_8422_2325, |h, &b| (h ^ u64::from(b)).wrapping_mul(0x0100_0000_01B3))
}

/// A value in `lo..=hi` from `r`.
fn within((lo, hi): (u32, u32), r: u64) -> u32 {
    lo + (r % (u64::from(hi - lo) + 1)) as u32
}

impl NameShape {
    /// The shape most of `names` have, and how many fell in each class.
    ///
    /// A name is 7-part or 5-part when Picard reads a position from it, and
    /// other otherwise. The shape is the class most names fall in; a tie, or no
    /// names, gives `Other`. Its middle is the most common one among names of
    /// that class (a tie goes to the smallest), and its tiles and x and y
    /// ranges are those of the names with that middle.
    pub fn learn<'a>(names: impl IntoIterator<Item = &'a [u8]>) -> (NameShape, NameClasses) {
        let mut classes = NameClasses::default();
        let (mut seven, mut five) = (Vec::new(), Vec::new());
        for name in names {
            match positioned(name) {
                Some(p) if p.fields == 7 => {
                    classes.seven += 1;
                    seven.push(p);
                }
                Some(p) => {
                    classes.five += 1;
                    five.push(p);
                }
                None => classes.other += 1,
            }
        }
        let class = if classes.seven > classes.five && classes.seven > classes.other {
            seven
        } else if classes.five > classes.seven && classes.five > classes.other {
            five
        } else {
            return (NameShape::Other, classes);
        };

        let mut counts: std::collections::BTreeMap<&[u8], usize> = std::collections::BTreeMap::new();
        for p in &class {
            *counts.entry(p.middle).or_default() += 1;
        }
        // In ascending order, so only a strictly larger count replaces the
        // best: a tie keeps the smallest middle.
        let mut middle: &[u8] = &[];
        let mut best = 0;
        for (&m, &n) in &counts {
            if n > best {
                (middle, best) = (m, n);
            }
        }
        let with_middle: Vec<&Positioned> = class.iter().filter(|p| p.middle == middle).collect();
        let tiles: std::collections::BTreeSet<u32> = with_middle.iter().map(|p| p.tile).collect();
        let range = |v: &dyn Fn(&Positioned) -> u32| {
            with_middle.iter().fold((u32::MAX, 0), |(lo, hi), p| (lo.min(v(p)), hi.max(v(p))))
        };
        let shape = NameShape::Illumina {
            middle: String::from_utf8_lossy(middle).into_owned(),
            tiles: tiles.into_iter().collect(),
            x: range(&|p| p.x),
            y: range(&|p| p.y),
        };
        (shape, classes)
    }

    /// The name spike writes for its read `internal` (e.g. `ev0001_hap_000123`):
    /// `SPIKE_<internal>`, then for an Illumina shape `:<middle>:TILE:X:Y`.
    ///
    /// The tile, x and y come from a hash of `internal`, never from the run's
    /// random stream, so naming changes nothing else about a run.
    pub fn name(&self, internal: &str) -> String {
        match self {
            NameShape::Other => format!("{MARK}{internal}"),
            NameShape::Illumina { middle, tiles, x, y } => {
                let mut state = fnv1a64(internal.as_bytes());
                let tile = tiles[(splitmix64(&mut state) % tiles.len() as u64) as usize];
                let px = within(*x, splitmix64(&mut state));
                let py = within(*y, splitmix64(&mut state));
                format!("{MARK}{internal}:{middle}:{tile}:{px}:{py}")
            }
        }
    }

    /// The shape in a few words, for the log.
    pub fn describe(&self) -> String {
        match self {
            NameShape::Other => "no flowcell position Picard can read".to_string(),
            NameShape::Illumina { middle, tiles, x, y } => format!(
                "MACHINE:{middle}:TILE:X:Y, {} tiles, x {}-{}, y {}-{}",
                tiles.len(),
                x.0,
                x.1,
                y.0,
                y.1
            ),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn learn(names: &[&str]) -> (NameShape, NameClasses) {
        NameShape::learn(names.iter().map(|n| n.as_bytes()))
    }

    /// Picard 3.3.0's optimized rule: 5 or 7 `:` fields, the last three
    /// numbers (tile, x, y).
    fn picard_position(name: &str) -> Option<(u32, u32, u32)> {
        let f: Vec<&str> = name.split(':').collect();
        let n = f.len();
        if n != 5 && n != 7 {
            return None;
        }
        Some((f[n - 3].parse().ok()?, f[n - 2].parse().ok()?, f[n - 1].parse().ok()?))
    }

    fn illumina(middle: &str, tiles: &[u32], x: (u32, u32), y: (u32, u32)) -> NameShape {
        NameShape::Illumina { middle: middle.to_string(), tiles: tiles.to_vec(), x, y }
    }

    #[test]
    fn test_a_five_or_seven_field_name_with_a_numeric_tail_has_a_position() {
        assert_eq!(
            positioned(b"A00744:46:HV3C3DSXX:2:1221:8775:9361"),
            Some(Positioned { fields: 7, middle: b"46:HV3C3DSXX:2", tile: 1221, x: 8775, y: 9361 })
        );
        assert_eq!(
            positioned(b"HWUSI-EAS100R:6:73:941:1973"),
            Some(Positioned { fields: 5, middle: b"6", tile: 73, x: 941, y: 1973 })
        );
        assert_eq!(
            positioned(b"M:1:FC:2:1101:123456789:5").map(|p| p.x),
            Some(123_456_789),
            "nine digits fit Picard's int"
        );
    }

    #[test]
    fn test_any_other_name_has_no_position() {
        for name in [
            "a:b:c:1:2:3",                  // 6 fields
            "a:b:c:d:e:f:1:2",              // 8 fields
            "M:1:FC:2:1101:12a:5",          // a letter in x
            "M:1:FC:2:1101::5",             // an empty x
            "M:1:FC:2:1101:+5:5",           // a sign
            "M:1:FC:2:1101:1234567890:5",   // ten digits
            "SRR062634.1",
            "chrA_pair0",
            "",
        ] {
            assert_eq!(positioned(name.as_bytes()), None, "{name}");
        }
    }

    #[test]
    fn test_the_shape_is_the_class_most_names_fall_in() {
        let (shape, classes) =
            learn(&["M:1:FC:2:1101:10:20", "M:1:FC:2:1101:11:21", "M:1:FC:2:1102:12:22", "r1", "r2"]);
        assert_eq!(classes, NameClasses { seven: 3, five: 0, other: 2 });
        assert_eq!(shape, illumina("1:FC:2", &[1101, 1102], (10, 12), (20, 22)));

        let (shape, classes) = learn(&["M:3:1101:10:20", "M:3:1101:12:22", "M:1:FC:2:1101:11:21", "r1"]);
        assert_eq!(classes, NameClasses { seven: 1, five: 2, other: 1 });
        assert_eq!(shape, illumina("3", &[1101], (10, 12), (20, 22)));
    }

    #[test]
    fn test_a_tie_or_no_names_gives_other() {
        let (shape, classes) = learn(&["M:1:FC:2:1101:10:20", "M:1:FC:2:1101:11:21", "r1", "r2"]);
        assert_eq!(classes, NameClasses { seven: 2, five: 0, other: 2 });
        assert_eq!(shape, NameShape::Other, "seven ties other");
        let (shape, _) = learn(&["M:1:FC:2:1101:10:20", "M:1:FC:2:1101:11:21", "M:3:1101:10:20", "M:3:1101:12:22"]);
        assert_eq!(shape, NameShape::Other, "seven ties five");
        assert_eq!(learn(&[]), (NameShape::Other, NameClasses::default()));
    }

    #[test]
    fn test_the_middle_is_the_most_common_one_and_a_tie_goes_to_the_smallest() {
        // The first name holds the minority middle.
        let (shape, _) = learn(&["M:1:FC:1:1101:5:5", "M:1:FC:2:1101:6:6", "M:1:FC:2:1102:7:7"]);
        assert_eq!(shape, illumina("1:FC:2", &[1101, 1102], (6, 7), (6, 7)));
        // A tie, the larger middle first.
        let (shape, _) = learn(&["M:1:FC:2:1101:6:6", "M:1:FC:1:1103:5:5"]);
        assert_eq!(shape, illumina("1:FC:1", &[1103], (5, 5), (5, 5)));
    }

    #[test]
    fn test_tiles_and_ranges_come_only_from_names_with_that_middle() {
        let (shape, _) = learn(&[
            "M:1:FC:2:1101:100:200",
            "M:1:FC:2:1102:300:50",
            "M:1:FC:2:1101:150:150",
            "M:1:FC:1:2201:99999:1", // another lane
        ]);
        assert_eq!(shape, illumina("1:FC:2", &[1101, 1102], (100, 300), (50, 200)));
    }

    #[test]
    fn test_an_other_shape_name_is_marked_and_has_no_position() {
        let name = NameShape::Other.name("ev0001_hap_000001");
        assert_eq!(name, "SPIKE_ev0001_hap_000001");
        assert_eq!(picard_position(&name), None);
    }

    #[test]
    fn test_an_illumina_shape_name_is_marked_and_picard_reads_a_position_in_range() {
        let tiles = [1101, 1102, 1103];
        let shape = illumina("46:FC:2", &tiles, (1000, 2000), (5000, 5100));
        let mut used = std::collections::BTreeSet::new();
        let (mut x_lo, mut x_hi) = (u32::MAX, 0);
        for i in 0..1000 {
            let internal = format!("ev0001_hap_{:06}", i);
            let name = shape.name(&internal);
            assert!(name.starts_with(&format!("SPIKE_{internal}:46:FC:2:")), "{name}");
            let (tile, x, y) = picard_position(&name).unwrap_or_else(|| panic!("no position: {name}"));
            assert!(tiles.contains(&tile), "{name}");
            assert!((1000..=2000).contains(&x) && (5000..=5100).contains(&y), "{name}");
            used.insert(tile);
            (x_lo, x_hi) = (x_lo.min(x), x_hi.max(x));
        }
        assert_eq!(used.len(), 3, "every tile is used");
        assert!(x_lo < 1100 && x_hi > 1900, "x spreads over its range: {x_lo}-{x_hi}");

        let five = illumina("6", &[73], (941, 941), (1973, 1973));
        assert_eq!(five.name("ev0002_dup_depth_000004"), "SPIKE_ev0002_dup_depth_000004:6:73:941:1973");
        assert_eq!(picard_position(&five.name("ev0002_dup_depth_000004")), Some((73, 941, 1973)));
    }

    #[test]
    fn test_the_same_internal_name_always_gets_the_same_name() {
        let shape = illumina("46:FC:2", &[1101, 1102, 1103], (1000, 2000), (5000, 5100));
        assert_eq!(shape.name("ev0001_hap_000007"), shape.clone().name("ev0001_hap_000007"));
        assert_ne!(shape.name("ev0001_hap_000007"), shape.name("ev0001_hap_000008"));
    }
}
