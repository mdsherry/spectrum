#[derive(Debug, PartialEq)]
pub struct RGB(u8, u8, u8);
#[derive(Debug, PartialEq)]
struct RGBf(f64, f64, f64);
#[derive(Debug, PartialEq)]
pub struct HSL(i16, f64, f64);
#[derive(Debug, PartialEq)]
pub struct HSV(i16, f64, f64);
#[derive(Debug, PartialEq)]
pub struct CMYK(f64, f64, f64, f64);
#[derive(Debug, PartialEq)]
pub struct YCbCr(u8, u8, u8);
#[derive(Debug, PartialEq)]
pub struct YCbCr709(u8, u8, u8);
#[derive(Debug, PartialEq)]
pub struct YCbCr2020(u8, u8, u8);

#[derive(Debug, PartialEq)]
pub struct YCbCrJPEG(u8, u8, u8);

macro_rules! chain_convert {
    ($init:ty; $($intermediate: ty),+; $dest: ty) => {
      impl From<$init> for $dest {
        fn from(x: $init) -> Self {
          $(
            let x: $intermediate = x.into();
          )*
          x.into()
        }
      }
    }
}


chain_convert!(HSL; RGBf; RGB);
chain_convert!(RGB; RGBf; HSL);

chain_convert!(HSV; HSL, RGBf; RGB);
chain_convert!(RGB; RGBf, HSL; HSV);

chain_convert!(CMYK; RGBf; RGB);
chain_convert!(RGB; RGBf; CMYK);

chain_convert!(CMYK; RGBf; HSL);
chain_convert!(CMYK; RGBf, HSL; HSV);

chain_convert!(HSL; RGBf; CMYK);
chain_convert!(HSV; HSL, RGBf; CMYK);

impl From<RGB> for RGBf {
    fn from(rgb: RGB) -> Self {
        let RGB(r, g, b) = rgb;
        let (r, g, b) = (r as f64 / 255., g as f64 / 255., b as f64 / 255.);
        RGBf(r, g, b)
    }
}

impl From<RGBf> for RGB {
    fn from(rgbf: RGBf) -> Self {
        let RGBf(r, g, b) = rgbf;
        let epsilon = 0.5;
        RGB((r * 255. + epsilon) as u8,
            (g * 255. + epsilon) as u8,
            (b * 255. + epsilon) as u8)
    }
}

impl From<HSL> for RGBf {
    fn from(hsl: HSL) -> Self {
        let HSL(h, s, l) = hsl;
        let c = (1. - (2. * l - 1.).abs()) * s;
        let h_prime = (h as f64) / 60.;
        let x = c * (1. - ((h_prime % 2.) - 1.).abs());
        let (r, g, b) = if h < 60 {
            (c, x, 0.)
        } else if h < 120 {
            (x, c, 0.)
        } else if h < 180 {
            (0., c, x)
        } else if h < 240 {
            (0., x, c)
        } else if h < 300 {
            (x, 0., c)
        } else {
            (c, 0., x)
        };
        let m = l - c / 2.;
        let (r, g, b) = (r + m, g + m, b + m);
        RGBf(r, g, b)
    }
}

impl From<RGBf> for HSL {
    fn from(rgbf: RGBf) -> Self {
        let RGBf(r, g, b) = rgbf;
        let cmax = r.max(g.max(b));
        let cmin = r.min(g.min(b));
        let delta = cmax - cmin;
        let h = if delta == 0. {
            // Black/grey/white
            0.
        } else if cmax == r {
            60. * (((g - b) / delta) % 6.)
        } else if cmax == g {
            60. * (((b - r) / delta) + 2.)
        } else {
            60. * (((r - g) / delta) + 4.)
        };
        let h = if h < 0. { h + 360. } else { h };
        let l = (cmax + cmin) / 2.;
        let s = if delta == 0. {
            0.
        } else {
            delta / (1. - (2. * l - 1.).abs())
        };
        HSL(h as i16, s, l)
    }
}


impl From<HSL> for HSV {
    fn from(hsl: HSL) -> Self {
        let HSL(h, s_hsl, l) = hsl;
        let v = (2. * l + s_hsl * (1. - (2. * l - 1.).abs())) / 2.;
        let s_hsv = 2. * (v - l) / v;
        HSV(h, s_hsv, v)
    }
}


impl From<HSV> for HSL {
    fn from(hsv: HSV) -> Self {
        let HSV(h, s_hsv, v) = hsv;
        let l = v * (2. * s_hsv) / 2.;
        let s_hsl = v * s_hsv / (1. - (2. * l - 1.).abs());
        HSL(h, s_hsl, l)
    }
}

impl From<CMYK> for RGBf {
    fn from(cmyk: CMYK) -> Self {
        let CMYK(c, m, y, k) = cmyk;
        let r = (1. - c) * (1. - k);
        let g = (1. - m) * (1. - k);
        let b = (1. - y) * (1. - k);
        RGBf(r, g, b)
    }
}

impl From<RGBf> for CMYK {
    fn from(rgbf: RGBf) -> Self {
        let RGBf(r, g, b) = rgbf;
        let k = 1. - r.max(g.max(b));
        if k == 1. {
            CMYK(0., 0., 0., 1.)
        } else {
            let c = (1. - r - k) / (1. - k);
            let m = (1. - g - k) / (1. - k);
            let y = (1. - b - k) / (1. - k);
            CMYK(c, m, y, k)
        }
    }
}

fn convert_rgb2yCrCb(k_r: f64, k_b: f64, r: f64, g: f64, b: f64) -> (u8, u8, u8) {
    let k_g = 1. - k_r - k_b;
    let y = r * k_r + g * k_g + b * k_b;
    let p_r = (r - y) / (2. * (1. - k_r));
    let p_b = (b - y) / (2. * (1. - k_b));
    let y = 16 + (y * 255.) as u8;
    let cr = 128 + (p_r * 255.) as u8;
    let cb = 128 + (p_b * 255.) as u8;
    (y, cb, cr)
}

impl From<RGBf> for YCbCr {
    fn from(rgbf: RGBf) -> Self {
        let RGBf(r, g, b) = rgbf;
        let (y, cb, cr) = convert_rgb2yCrCb(0.299, 0.114, r, g, b);
        YCbCr(y, cb, cr)
    }
}

impl From<RGBf> for YCbCr709 {
    fn from(rgbf: RGBf) -> Self {
        let RGBf(r, g, b) = rgbf;
        let (y, cb, cr) = convert_rgb2yCrCb(0.2126, 0.0722, r, g, b);
        YCbCr709(y, cb, cr)
    }
}

impl From<RGBf> for YCbCr2020 {
    fn from(rgbf: RGBf) -> Self {
        let RGBf(r, g, b) = rgbf;
        let (y, cb, cr) = convert_rgb2yCrCb(0.2627, 0.0593, r, g, b);
        YCbCr2020(y, cb, cr)
    }
}


impl From<RGB> for YCbCrJPEG {
    fn from(rgb: RGB) -> Self {
        let RGB(r, g, b) = rgb;
        let r = r as f64;
        let g = g as f64;
        let b = b as f64;
        let y = 0.299 * r + 0.587 * g + 0.114 * b;
        let cb = 128. + (0.168736 * r - 0.331264 * g + 0.5 * b);
        let cr = 128. + (0.5 * r - 0.418688 * g - 0.081312 * b);
        YCbCrJPEG(y as u8, cb as u8, cr as u8)
    }
}

#[cfg(test)]
mod tests {
    use HSL;
    use RGB;
    #[test]
    fn hsl_to_rgb_works() {
        // Black
        assert_eq!(RGB(0, 0, 0), HSL(0, 0., 0.).into());
        // White
        assert_eq!(RGB(255, 255, 255), HSL(0, 0., 1.).into());
        // Red
        assert_eq!(RGB(255, 0, 0), HSL(0, 1., 0.5).into());
        // Lime
        assert_eq!(RGB(0, 255, 0), HSL(120, 1., 0.5).into());
        // Blue
        assert_eq!(RGB(0, 0, 255), HSL(240, 1., 0.5).into());
        // Yellow
        assert_eq!(RGB(255, 255, 0), HSL(60, 1., 0.5).into());
        // Cyan
        assert_eq!(RGB(0, 255, 255), HSL(180, 1., 0.5).into());
        // Magenta
        assert_eq!(RGB(255, 0, 255), HSL(300, 1., 0.5).into());
        // Silver
        assert_eq!(RGB(192, 192, 192), HSL(0, 0., 0.753).into());
        // Gray
        assert_eq!(RGB(128, 128, 128), HSL(0, 0., 0.5).into());
        // Maroon
        assert_eq!(RGB(128, 0, 0), HSL(0, 1., 0.25).into());
        // Olive
        assert_eq!(RGB(128, 128, 0), HSL(60, 1., 0.25).into());
        // Green
        assert_eq!(RGB(0, 128, 0), HSL(120, 1., 0.25).into());
        // Purple
        assert_eq!(RGB(128, 0, 128), HSL(300, 1., 0.25).into());
        // Teal
        assert_eq!(RGB(0, 128, 128), HSL(180, 1., 0.25).into());
        // Navy
        assert_eq!(RGB(0, 0, 128), HSL(240, 1., 0.25).into());

        assert_eq!(RGB(151, 195, 183), HSL(163, 0.27, 0.68).into());
    }

    #[test]
    fn rgb_to_hsl_works() {
        // Black
        assert_eq!(HSL(0, 0., 0.), RGB(0, 0, 0).into());
        // White
        assert_eq!(HSL(0, 0., 1.), RGB(255, 255, 255).into());
        // Red
        assert_eq!(HSL(0, 1., 0.5), RGB(255, 0, 0).into());
        // Lime
        assert_eq!(HSL(120, 1., 0.5), RGB(0, 255, 0).into());
        // Blue
        assert_eq!(HSL(240, 1., 0.5), RGB(0, 0, 255).into());
        // Yellow
        assert_eq!(HSL(60, 1., 0.5), RGB(255, 255, 0).into());
        // Cyan
        assert_eq!(HSL(180, 1., 0.5), RGB(0, 255, 255).into());
        // Magenta
        assert_eq!(HSL(300, 1., 0.5), RGB(255, 0, 255).into());
        // Silver
        assert_eq!(HSL(0, 0., 0.7529411764705882), RGB(192, 192, 192).into());
        // Gray
        assert_eq!(HSL(0, 0., 0.5019607843137255), RGB(128, 128, 128).into());
        // Maroon
        assert_eq!(HSL(0, 1., 0.25098039215686274), RGB(128, 0, 0).into());
        // Olive
        assert_eq!(HSL(60, 1., 0.25098039215686274), RGB(128, 128, 0).into());
        // Green
        assert_eq!(HSL(120, 1., 0.25098039215686274), RGB(0, 128, 0).into());
        // Purple
        assert_eq!(HSL(300, 1., 0.25098039215686274), RGB(128, 0, 128).into());
        // Teal
        assert_eq!(HSL(180, 1., 0.25098039215686274), RGB(0, 128, 128).into());
        // Navy
        assert_eq!(HSL(240, 1., 0.25098039215686274), RGB(0, 0, 128).into());

        assert_eq!(HSL(163, 0.26829268292682923, 0.6784313725490196),
                   RGB(151, 195, 183).into());
    }
}
