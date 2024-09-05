//! This module contains the iterators returned by various functions on the caches of bessel function values.

// These iterators are simple wrappers of the underlying `HashMap` iterators that hide the implementation detail that we actually hash `u64`, and not `f64`.

pub use cached_jnu_ynu::{CachedJnuYnuIter, CachedJnuYnuArguments, CachedJnuYnuValues};
mod cached_jnu_ynu {
    use std::{collections::hash_map::Iter, iter::FusedIterator};
    /// An iterator created by the [`CachedJnuYnu::iter`](super::super::CachedJnuYnu::iter) function, see it for more information.
    #[derive(Debug, Clone)]
    pub struct CachedJnuYnuIter<'a>(
        core::iter::Map<std::collections::hash_map::Iter<'a, (u64, u64), (f64, f64)>, fn((&(u64, u64), &(f64, f64))) -> ((f64, f64), (f64, f64))>
    );

    impl<'a> CachedJnuYnuIter<'a> {
        pub(crate) fn new(cache_iter: Iter<'a, (u64, u64), (f64, f64)>) -> Self {
            Self(cache_iter.map(|(k, v)| ((f64::from_bits(k.0), f64::from_bits(k.1)), *v)))
        }
    }

    impl<'a> Iterator for CachedJnuYnuIter<'a> {
        type Item = ((f64, f64), (f64, f64));

        #[inline]
        fn next(&mut self) -> Option<Self::Item> {
            self.0.next()
        }

        #[inline]
        fn size_hint(&self) -> (usize, Option<usize>) {
            self.0.size_hint()
        }

        #[inline]
        fn nth(&mut self, n: usize) -> Option<Self::Item> {
            self.0.nth(n)
        }

        #[inline]
        fn last(self) -> Option<Self::Item> {
            self.0.last()
        }

        #[inline]
        fn count(self) -> usize {
            self.0.count()
        }
    }
    impl<'a> ExactSizeIterator for CachedJnuYnuIter<'a> {
        #[inline]
        fn len(&self) -> usize {
            self.0.len()
        }
    }
    impl<'a> FusedIterator for CachedJnuYnuIter<'a> {}

    /// An iterator created by the [`CachedJnuYnu::arguments`](super::super::CachedJnuYnu::arguments) function, see it for more information.
    #[derive(Debug, Clone)]
    pub struct CachedJnuYnuArguments<'a>(core::iter::Map<std::collections::hash_map::Keys<'a, (u64, u64), (f64, f64)>, fn(&(u64, u64)) -> (f64, f64)>);

    impl<'a> CachedJnuYnuArguments<'a> {
        pub(crate) fn new(cache_iter: std::collections::hash_map::Keys<'a, (u64, u64), (f64, f64)>) -> Self {
            Self(
                cache_iter.map(|k| (f64::from_bits(k.0), f64::from_bits(k.1)))
            )
        }
    }

    impl<'a> Iterator for CachedJnuYnuArguments<'a> {
        type Item = (f64, f64);

        #[inline]
        fn next(&mut self) -> Option<Self::Item> {
            self.0.next()
        }

        #[inline]
        fn size_hint(&self) -> (usize, Option<usize>) {
            self.0.size_hint()
        }

        #[inline]
        fn nth(&mut self, n: usize) -> Option<Self::Item> {
            self.0.nth(n)
        }

        #[inline]
        fn last(self) -> Option<Self::Item> {
            self.0.last()
        }

        #[inline]
        fn count(self) -> usize {
            self.0.count()
        }
    }
    impl<'a> ExactSizeIterator for CachedJnuYnuArguments<'a> {
        #[inline]
        fn len(&self) -> usize {
            self.0.len()
        }
    }
    impl<'a> FusedIterator for CachedJnuYnuArguments<'a> {}

    /// An iterator created by the [`CachedJnuYnu::values`](super::super::CachedJnuYnu::values) function, see it for more information.
    #[derive(Debug, Clone)]
    pub struct CachedJnuYnuValues<'a>(core::iter::Copied<std::collections::hash_map::Values<'a, (u64, u64), (f64, f64)>>);

    impl<'a> CachedJnuYnuValues<'a> {
        pub(crate) fn new(cache_iter: std::collections::hash_map::Values<'a, (u64, u64), (f64, f64)>) -> Self {
            Self(
                cache_iter.copied()
            )
        }
    }

    impl<'a> Iterator for CachedJnuYnuValues<'a> {
        type Item = (f64, f64);

        #[inline]
        fn next(&mut self) -> Option<Self::Item> {
            self.0.next()
        }

        #[inline]
        fn size_hint(&self) -> (usize, Option<usize>) {
            self.0.size_hint()
        }

        #[inline]
        fn nth(&mut self, n: usize) -> Option<Self::Item> {
            self.0.nth(n)
        }

        #[inline]
        fn last(self) -> Option<Self::Item> {
            self.0.last()
        }

        #[inline]
        fn count(self) -> usize {
            self.0.count()
        }
    }
    impl<'a> ExactSizeIterator for CachedJnuYnuValues<'a> {
        #[inline]
        fn len(&self) -> usize {
            self.0.len()
        }
    }
    impl<'a> FusedIterator for CachedJnuYnuValues<'a> {}
}

pub use cached_inu_knu::{CachedInuKnuIter, CachedInuKnuArguments, CachedInuKnuValues};
mod cached_inu_knu {
    use std::{collections::hash_map::Iter, iter::FusedIterator};
    /// An iterator created by the [`CachedInuKnu::iter`](super::super::CachedInuKnu::iter) function, see it for more information.
    #[derive(Debug, Clone)]
    pub struct CachedInuKnuIter<'a>(
        core::iter::Map<std::collections::hash_map::Iter<'a, (u64, u64), (f64, f64)>, fn((&(u64, u64), &(f64, f64))) -> ((f64, f64), (f64, f64))>
    );

    impl<'a> CachedInuKnuIter<'a> {
        pub(crate) fn new(cache_iter: Iter<'a, (u64, u64), (f64, f64)>) -> Self {
            Self(cache_iter.map(|(k, v)| ((f64::from_bits(k.0), f64::from_bits(k.1)), *v)))
        }
    }

    impl<'a> Iterator for CachedInuKnuIter<'a> {
        type Item = ((f64, f64), (f64, f64));

        #[inline]
        fn next(&mut self) -> Option<Self::Item> {
            self.0.next()
        }

        #[inline]
        fn size_hint(&self) -> (usize, Option<usize>) {
            self.0.size_hint()
        }

        #[inline]
        fn nth(&mut self, n: usize) -> Option<Self::Item> {
            self.0.nth(n)
        }

        #[inline]
        fn last(self) -> Option<Self::Item> {
            self.0.last()
        }

        #[inline]
        fn count(self) -> usize {
            self.0.count()
        }
    }
    impl<'a> ExactSizeIterator for CachedInuKnuIter<'a> {
        #[inline]
        fn len(&self) -> usize {
            self.0.len()
        }
    }
    impl<'a> FusedIterator for CachedInuKnuIter<'a> {}

    /// An iterator created by the [`CachedInuKnu::arguments`](super::super::CachedInuKnu::arguments) function, see it for more information.
    #[derive(Debug, Clone)]
    pub struct CachedInuKnuArguments<'a>(core::iter::Map<std::collections::hash_map::Keys<'a, (u64, u64), (f64, f64)>, fn(&(u64, u64)) -> (f64, f64)>);

    impl<'a> CachedInuKnuArguments<'a> {
        pub(crate) fn new(cache_iter: std::collections::hash_map::Keys<'a, (u64, u64), (f64, f64)>) -> Self {
            Self(
                cache_iter.map(|k| (f64::from_bits(k.0), f64::from_bits(k.1)))
            )
        }
    }

    impl<'a> Iterator for CachedInuKnuArguments<'a> {
        type Item = (f64, f64);

        #[inline]
        fn next(&mut self) -> Option<Self::Item> {
            self.0.next()
        }

        #[inline]
        fn size_hint(&self) -> (usize, Option<usize>) {
            self.0.size_hint()
        }

        #[inline]
        fn nth(&mut self, n: usize) -> Option<Self::Item> {
            self.0.nth(n)
        }

        #[inline]
        fn last(self) -> Option<Self::Item> {
            self.0.last()
        }

        #[inline]
        fn count(self) -> usize {
            self.0.count()
        }
    }
    impl<'a> ExactSizeIterator for CachedInuKnuArguments<'a> {
        #[inline]
        fn len(&self) -> usize {
            self.0.len()
        }
    }
    impl<'a> FusedIterator for CachedInuKnuArguments<'a> {}

    /// An iterator created by the [`CachedInuKnu::values`](super::super::CachedInuKnu::values) function, see it for more information.
    #[derive(Debug, Clone)]
    pub struct CachedInuKnuValues<'a>(core::iter::Copied<std::collections::hash_map::Values<'a, (u64, u64), (f64, f64)>>);

    impl<'a> CachedInuKnuValues<'a> {
        pub(crate) fn new(cache_iter: std::collections::hash_map::Values<'a, (u64, u64), (f64, f64)>) -> Self {
            Self(
                cache_iter.copied()
            )
        }
    }

    impl<'a> Iterator for CachedInuKnuValues<'a> {
        type Item = (f64, f64);

        #[inline]
        fn next(&mut self) -> Option<Self::Item> {
            self.0.next()
        }

        #[inline]
        fn size_hint(&self) -> (usize, Option<usize>) {
            self.0.size_hint()
        }

        #[inline]
        fn nth(&mut self, n: usize) -> Option<Self::Item> {
            self.0.nth(n)
        }

        #[inline]
        fn last(self) -> Option<Self::Item> {
            self.0.last()
        }

        #[inline]
        fn count(self) -> usize {
            self.0.count()
        }
    }
    impl<'a> ExactSizeIterator for CachedInuKnuValues<'a> {
        #[inline]
        fn len(&self) -> usize {
            self.0.len()
        }
    }
    impl<'a> FusedIterator for CachedInuKnuValues<'a> {}
}

pub use cached_besselik::{CachedBesselIKIter, CachedBesselIKArguments, CachedBesselIKValues};
mod cached_besselik {
    use std::{collections::hash_map::Iter, iter::FusedIterator};
    /// An iterator created by the [`CachedBesselIK::iter`](super::super::CachedBesselIK::iter) function, see it for more information.
    #[derive(Debug, Clone)]
    pub struct CachedBesselIKIter<'a>(
        core::iter::Map<std::collections::hash_map::Iter<'a, (u64, u64), (f64, f64, f64, f64)>, fn((&(u64, u64), &(f64, f64, f64, f64))) -> ((f64, f64), (f64, f64, f64, f64))>
    );

    impl<'a> CachedBesselIKIter<'a> {
        pub(crate) fn new(cache_iter: Iter<'a, (u64, u64), (f64, f64, f64, f64)>) -> Self {
            Self(cache_iter.map(|(k, v)| ((f64::from_bits(k.0), f64::from_bits(k.1)), *v)))
        }
    }

    impl<'a> Iterator for CachedBesselIKIter<'a> {
        type Item = ((f64, f64), (f64, f64, f64, f64));

        #[inline]
        fn next(&mut self) -> Option<Self::Item> {
            self.0.next()
        }

        #[inline]
        fn size_hint(&self) -> (usize, Option<usize>) {
            self.0.size_hint()
        }

        #[inline]
        fn nth(&mut self, n: usize) -> Option<Self::Item> {
            self.0.nth(n)
        }

        #[inline]
        fn last(self) -> Option<Self::Item> {
            self.0.last()
        }

        #[inline]
        fn count(self) -> usize {
            self.0.count()
        }
    }
    impl<'a> ExactSizeIterator for CachedBesselIKIter<'a> {
        #[inline]
        fn len(&self) -> usize {
            self.0.len()
        }
    }
    impl<'a> FusedIterator for CachedBesselIKIter<'a> {}

    /// An iterator created by the [`CachedBesselIK::arguments`](super::super::CachedBesselIK::arguments) function, see it for more information.
    #[derive(Debug, Clone)]
    pub struct CachedBesselIKArguments<'a>(core::iter::Map<std::collections::hash_map::Keys<'a, (u64, u64), (f64, f64, f64, f64)>, fn(&(u64, u64)) -> (f64, f64)>);

    impl<'a> CachedBesselIKArguments<'a> {
        pub(crate) fn new(cache_iter: std::collections::hash_map::Keys<'a, (u64, u64), (f64, f64, f64, f64)>) -> Self {
            Self(
                cache_iter.map(|k| (f64::from_bits(k.0), f64::from_bits(k.1)))
            )
        }
    }

    impl<'a> Iterator for CachedBesselIKArguments<'a> {
        type Item = (f64, f64);

        #[inline]
        fn next(&mut self) -> Option<Self::Item> {
            self.0.next()
        }

        #[inline]
        fn size_hint(&self) -> (usize, Option<usize>) {
            self.0.size_hint()
        }

        #[inline]
        fn nth(&mut self, n: usize) -> Option<Self::Item> {
            self.0.nth(n)
        }

        #[inline]
        fn last(self) -> Option<Self::Item> {
            self.0.last()
        }

        #[inline]
        fn count(self) -> usize {
            self.0.count()
        }
    }
    impl<'a> ExactSizeIterator for CachedBesselIKArguments<'a> {
        #[inline]
        fn len(&self) -> usize {
            self.0.len()
        }
    }
    impl<'a> FusedIterator for CachedBesselIKArguments<'a> {}

    /// An iterator created by the [`CachedBesselIK::values`](super::super::CachedBesselIK::values) function, see it for more information.
    #[derive(Debug, Clone)]
    pub struct CachedBesselIKValues<'a>(core::iter::Copied<std::collections::hash_map::Values<'a, (u64, u64), (f64, f64, f64, f64)>>);

    impl<'a> CachedBesselIKValues<'a> {
        pub(crate) fn new(cache_iter: std::collections::hash_map::Values<'a, (u64, u64), (f64, f64, f64, f64)>) -> Self {
            Self(
                cache_iter.copied()
            )
        }
    }

    impl<'a> Iterator for CachedBesselIKValues<'a> {
        type Item = (f64, f64, f64, f64);

        #[inline]
        fn next(&mut self) -> Option<Self::Item> {
            self.0.next()
        }

        #[inline]
        fn size_hint(&self) -> (usize, Option<usize>) {
            self.0.size_hint()
        }

        #[inline]
        fn nth(&mut self, n: usize) -> Option<Self::Item> {
            self.0.nth(n)
        }

        #[inline]
        fn last(self) -> Option<Self::Item> {
            self.0.last()
        }

        #[inline]
        fn count(self) -> usize {
            self.0.count()
        }
    }
    impl<'a> ExactSizeIterator for CachedBesselIKValues<'a> {
        #[inline]
        fn len(&self) -> usize {
            self.0.len()
        }
    }
    impl<'a> FusedIterator for CachedBesselIKValues<'a> {}
}

pub use cached_besseljy::{CachedBesselJYIter, CachedBesselJYArguments, CachedBesselJYValues};
mod cached_besseljy {
    use std::{collections::hash_map::Iter, iter::FusedIterator};
    /// An iterator created by the [`CachedBesselJY::iter`](super::super::CachedBesselJY::iter) function, see it for more information.
    #[derive(Debug, Clone)]
    pub struct CachedBesselJYIter<'a>(
        core::iter::Map<std::collections::hash_map::Iter<'a, (u64, u64), (f64, f64, f64, f64)>, fn((&(u64, u64), &(f64, f64, f64, f64))) -> ((f64, f64), (f64, f64, f64, f64))>
    );

    impl<'a> CachedBesselJYIter<'a> {
        pub(crate) fn new(cache_iter: Iter<'a, (u64, u64), (f64, f64, f64, f64)>) -> Self {
            Self(cache_iter.map(|(k, v)| ((f64::from_bits(k.0), f64::from_bits(k.1)), *v)))
        }
    }

    impl<'a> Iterator for CachedBesselJYIter<'a> {
        type Item = ((f64, f64), (f64, f64, f64, f64));

        #[inline]
        fn next(&mut self) -> Option<Self::Item> {
            self.0.next()
        }

        #[inline]
        fn size_hint(&self) -> (usize, Option<usize>) {
            self.0.size_hint()
        }

        #[inline]
        fn nth(&mut self, n: usize) -> Option<Self::Item> {
            self.0.nth(n)
        }

        #[inline]
        fn last(self) -> Option<Self::Item> {
            self.0.last()
        }

        #[inline]
        fn count(self) -> usize {
            self.0.count()
        }
    }
    impl<'a> ExactSizeIterator for CachedBesselJYIter<'a> {
        #[inline]
        fn len(&self) -> usize {
            self.0.len()
        }
    }
    impl<'a> FusedIterator for CachedBesselJYIter<'a> {}

    /// An iterator created by the [`CachedBesselJY::arguments`](super::super::CachedBesselJY::arguments) function, see it for more information.
    #[derive(Debug, Clone)]
    pub struct CachedBesselJYArguments<'a>(core::iter::Map<std::collections::hash_map::Keys<'a, (u64, u64), (f64, f64, f64, f64)>, fn(&(u64, u64)) -> (f64, f64)>);

    impl<'a> CachedBesselJYArguments<'a> {
        pub(crate) fn new(cache_iter: std::collections::hash_map::Keys<'a, (u64, u64), (f64, f64, f64, f64)>) -> Self {
            Self(
                cache_iter.map(|k| (f64::from_bits(k.0), f64::from_bits(k.1)))
            )
        }
    }

    impl<'a> Iterator for CachedBesselJYArguments<'a> {
        type Item = (f64, f64);

        #[inline]
        fn next(&mut self) -> Option<Self::Item> {
            self.0.next()
        }

        #[inline]
        fn size_hint(&self) -> (usize, Option<usize>) {
            self.0.size_hint()
        }

        #[inline]
        fn nth(&mut self, n: usize) -> Option<Self::Item> {
            self.0.nth(n)
        }

        #[inline]
        fn last(self) -> Option<Self::Item> {
            self.0.last()
        }

        #[inline]
        fn count(self) -> usize {
            self.0.count()
        }
    }
    impl<'a> ExactSizeIterator for CachedBesselJYArguments<'a> {
        #[inline]
        fn len(&self) -> usize {
            self.0.len()
        }
    }
    impl<'a> FusedIterator for CachedBesselJYArguments<'a> {}

    /// An iterator created by the [`CachedBesselJY::values`](super::super::CachedBesselJY::values) function, see it for more information.
    #[derive(Debug, Clone)]
    pub struct CachedBesselJYValues<'a>(core::iter::Copied<std::collections::hash_map::Values<'a, (u64, u64), (f64, f64, f64, f64)>>);

    impl<'a> CachedBesselJYValues<'a> {
        pub(crate) fn new(cache_iter: std::collections::hash_map::Values<'a, (u64, u64), (f64, f64, f64, f64)>) -> Self {
            Self(
                cache_iter.copied()
            )
        }
    }

    impl<'a> Iterator for CachedBesselJYValues<'a> {
        type Item = (f64, f64, f64, f64);

        #[inline]
        fn next(&mut self) -> Option<Self::Item> {
            self.0.next()
        }

        #[inline]
        fn size_hint(&self) -> (usize, Option<usize>) {
            self.0.size_hint()
        }

        #[inline]
        fn nth(&mut self, n: usize) -> Option<Self::Item> {
            self.0.nth(n)
        }

        #[inline]
        fn last(self) -> Option<Self::Item> {
            self.0.last()
        }

        #[inline]
        fn count(self) -> usize {
            self.0.count()
        }
    }
    impl<'a> ExactSizeIterator for CachedBesselJYValues<'a> {
        #[inline]
        fn len(&self) -> usize {
            self.0.len()
        }
    }
    impl<'a> FusedIterator for CachedBesselJYValues<'a> {}
}