pub trait Ensure<T>: Sized {
    /// Predicate check where the caller provides a "container hint" that determines output type
    /// Pass `None` to get `Option<T>`; pass `Err(e)` to get `Result<T, E>`
    fn ensure<F, R>(self, predicate: F, on_fail: R) -> R
    where
        F: FnOnce(&Self) -> bool,
        R: EnsureFrom<Self>;

    fn ensure_map<R, F, ME>(self, predicate: F, map: ME) -> R
    where
        F: FnOnce(&T) -> bool,
        ME: FnOnce(&T) -> R,
        R: EnsureFrom<Self>;

    fn ensure_map_err<E, F, ME>(self, predicate: F, map_err: ME) -> Result<T, E>
    where
        F: FnOnce(&T) -> bool,
        ME: FnOnce(&T) -> E;
}

pub trait EnsureFrom<T>: Sized {
    fn from_ok(val: T) -> Self;
    fn from_fail(self) -> Self;
}

impl<T> EnsureFrom<T> for Option<T> {
    #[inline]
    fn from_ok(val: T) -> Self {
        Some(val)
    }

    #[inline]
    fn from_fail(self) -> Self {
        self
    }
}

impl<T, E> EnsureFrom<T> for Result<T, E> {
    #[inline]
    fn from_ok(val: T) -> Self {
        Ok(val)
    }

    #[inline]
    fn from_fail(self) -> Self {
        self
    }
}

impl<T> Ensure<T> for T {
    #[inline]
    fn ensure<F, R>(self, predicate: F, on_fail: R) -> R
    where
        F: FnOnce(&Self) -> bool,
        R: EnsureFrom<Self>,
    {
        if predicate(&self) {
            R::from_ok(self)
        } else {
            on_fail.from_fail()
        }
    }

    #[inline]
    fn ensure_map<R, F, ME>(self, predicate: F, map: ME) -> R
    where
        F: FnOnce(&T) -> bool,
        ME: FnOnce(&T) -> R,
        R: EnsureFrom<Self>,
    {
        if predicate(&self) {
            R::from_ok(self)
        } else {
            map(&self).from_fail()
        }
    }

    #[inline]
    fn ensure_map_err<E, F, ME>(self, predicate: F, map_err: ME) -> Result<T, E>
    where
        F: FnOnce(&T) -> bool,
        ME: FnOnce(&T) -> E,
    {
        if predicate(&self) {
            Ok(self)
        } else {
            Err(map_err(&self))
        }
    }
}
